#!/usr/bin/env python
"""Verify the per-cell GeoJSON against the cells parquet and the images, before REDSEA.

Every failure mode here is silent. REDSEA rasterizes the polygons into an int32 instance mask,
sums each channel per cell, and joins the result back on the QuPath detection UUID. If the ids
don't match it produces a smaller-than-expected table; if the coordinates are in microns or a
downsampled space, it produces a *plausible* table of wrong numbers. Neither raises.

Checks, in the order they would bite:
  1. one feature per cell, counts matching the parquet per donor
  2. feature ids are UUIDs, unique, and set-equal to the parquet's object_id
  3. geometry is Polygon/MultiPolygon and lies inside the image's FULL-RESOLUTION pixel bounds
  4. coordinates are pixels, not microns -- cross-checked against the parquet's centroid_x/y_um
     via PIXEL_SIZE_UM, which is the check that catches a downsampled or unit-converted export

Usage:
    conda activate scvi-env
    python verify_geojson.py                 # all donors present
    python verify_geojson.py --donor 112
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import tifffile

from codex_panel import DONORS, normalize_donor_id

GEOJSON_DIR = Path("/data/codex_phenocycler/redsea_scratch/geojson")
CELLS_DIR = Path("/data/codex_phenocycler/cells")
IMAGES_DIR = Path(
    "/home/smith6jt/panc_CODEX/Islet-Explorer/app/shiny_app/www/local_images"
)
PIXEL_SIZE_UM = 0.3774          # repo-wide constant (00_globals.R); OME PhysicalSizeX agrees
UUID_LEN = 36


def image_path(donor: str) -> Path | None:
    """local_images files are '<status>_<zero-padded case>.ome.tiff'; donor ids are unpadded."""
    for p in IMAGES_DIR.glob("*.ome.tiff"):
        if normalize_donor_id(p.name.replace(".ome.tiff", "").split("_")[-1]) == donor:
            return p
    return None


def check(donor: str, con: duckdb.DuckDBPyConnection, sample: int) -> list[str]:
    errs: list[str] = []
    gj = GEOJSON_DIR / f"cells__{donor.zfill(4) if len(donor) < 4 else donor}.geojson"
    if not gj.exists():                       # image name is zero-padded, donor id is not
        alt = GEOJSON_DIR / f"cells__{donor}.geojson"
        gj = alt if alt.exists() else gj
    if not gj.exists():
        return [f"missing GeoJSON: {gj}"]

    feats = json.load(gj.open())["features"]
    part = CELLS_DIR / f"donor_id={donor}"
    pq = con.execute(
        f"SELECT object_id, centroid_x_um, centroid_y_um FROM read_parquet('{part}/*.parquet')"
    ).fetchdf()

    # 1. counts
    if len(feats) != len(pq):
        errs.append(f"count mismatch: {len(feats):,} features vs {len(pq):,} parquet rows")

    # 2. ids
    gids = [f.get("id") for f in feats]
    if any(g is None for g in gids):
        errs.append("some features have no `id` — REDSEA cannot join them")
    gset = set(map(str, gids))
    if len(gset) != len(gids):
        errs.append(f"duplicate feature ids: {len(gids) - len(gset):,}")
    bad_uuid = [g for g in list(gset)[:1000] if len(str(g)) != UUID_LEN]
    if bad_uuid:
        errs.append(f"ids are not UUIDs, e.g. {bad_uuid[:3]}")
    pset = set(pq.object_id.astype(str))
    only_gj, only_pq = gset - pset, pset - gset
    if only_gj or only_pq:
        errs.append(
            f"id sets differ: {len(only_gj):,} only in GeoJSON, {len(only_pq):,} only in parquet"
        )

    # 3-4. geometry: pixel bounds and the micron cross-check
    img = image_path(donor)
    if img is None:
        errs.append(f"no local_images file for donor {donor}")
    else:
        shape = tifffile.TiffFile(img).series[0].shape       # (C, Y, X)
        H, W = shape[-2], shape[-1]
        rng = np.random.default_rng(0)
        idx = rng.choice(len(feats), size=min(sample, len(feats)), replace=False)
        cent, ids, bad_geom = [], [], 0
        for i in idx:
            f = feats[int(i)]
            g = f.get("geometry") or {}
            if g.get("type") not in ("Polygon", "MultiPolygon"):
                bad_geom += 1
                continue
            cs = g["coordinates"]
            pts = np.array(cs[0] if g["type"] == "Polygon" else cs[0][0], dtype=float)
            cent.append(pts.mean(axis=0))
            ids.append(str(f.get("id")))
        if bad_geom:
            errs.append(f"{bad_geom} sampled features are not Polygon/MultiPolygon")
        if cent:
            cent = np.array(cent)
            oob = ((cent[:, 0] < 0) | (cent[:, 0] > W) | (cent[:, 1] < 0) | (cent[:, 1] > H)).sum()
            if oob:
                errs.append(
                    f"{oob}/{len(cent)} sampled centroids fall outside the full-res image "
                    f"({W}x{H}) — export is not in full-resolution pixel space"
                )
            # The decisive unit check: parquet centroids are microns, GeoJSON must be pixels.
            m = pq.set_index(pq.object_id.astype(str)).reindex(ids)
            ok = m.centroid_x_um.notna().to_numpy()
            if ok.sum() >= 10:
                px = np.c_[m.centroid_x_um.to_numpy()[ok], m.centroid_y_um.to_numpy()[ok]] / PIXEL_SIZE_UM
                d = np.linalg.norm(cent[ok] - px, axis=1)
                med = float(np.median(d))
                if med > 5.0:
                    ratio = float(np.median(np.linalg.norm(cent[ok], axis=1) /
                                            np.maximum(np.linalg.norm(px, axis=1), 1e-9)))
                    errs.append(
                        f"polygon centroids disagree with parquet centroids by {med:.1f} px "
                        f"(median); implied scale {ratio:.4f}. 1.0 = correct; ~{PIXEL_SIZE_UM} "
                        f"= exported in microns; 0.5/0.25 = downsampled"
                    )
    return errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", action="append", help="repeatable; default = every donor present")
    ap.add_argument("--sample", type=int, default=2000, help="polygons sampled for geometry checks")
    args = ap.parse_args()

    donors = args.donor or [d for d in DONORS
                            if (GEOJSON_DIR / f"cells__{d.zfill(4) if len(d) < 4 else d}.geojson").exists()
                            or (GEOJSON_DIR / f"cells__{d}.geojson").exists()]
    if not donors:
        raise SystemExit(f"no GeoJSON found in {GEOJSON_DIR}; run export_cells_geojson.sh first")

    con = duckdb.connect()
    failed = {}
    for d in donors:
        errs = check(d, con, args.sample)
        print(f"{d:>6}  {'OK' if not errs else 'FAIL'}")
        for e in errs:
            print(f"        - {e}")
        if errs:
            failed[d] = errs

    print(f"\n{len(donors) - len(failed)}/{len(donors)} donors verified")
    if len(donors) < len(DONORS):
        print(f"NOT YET EXPORTED: {sorted(set(DONORS) - set(donors), key=int)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
