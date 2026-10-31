#!/usr/bin/env python3
"""Measure union-clipped areas from existing v10 structure annotations.

Overlapping same-class annotation fragments contribute physical area only once.
Uses every saved Islet/EndocrineCluster core, including unpaired geometric
objects. No source geometry, marker gate, image, or app data is changed.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import hashlib
import json
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import shapely
from shapely.affinity import scale
from shapely.geometry import shape
from shapely import STRtree

ROOT = Path(__file__).resolve().parents[2]
APP = ROOT / "data/development/uniform-20donor-20261022-v10"
OUT = ROOT / "output/pancfloyd_20donor_20261026/slides1_4/data"
CENSUS = ROOT / "output/h5_h10_original_plan_20261026/h8/full_object_census.parquet"
OLD = ROOT / "tmp/ufrc-release-20261022/artifacts/release-d55680a/data/summary/targets.parquet"
KINDS = ["Capillary", "Nerve", "Lymphatic"]
PX = .3774


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8*1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def donor(case, expected):
    path = APP / "gson" / f"{case:04d}.geojson.gz"
    with gzip.open(path, "rt") as f:
        features = json.load(f)["features"]
    cores, expanded, tissue = {}, {}, []
    parts = {kind: [] for kind in KINDS}
    invalid = 0
    for feature in features:
        p = feature.get("properties") or {}
        cl = (p.get("classification") or {}).get("name")
        if cl not in ["Islet", "IsletExpanded", "Tissue"] + KINDS:
            continue
        g = shape(feature["geometry"])
        if not g.is_valid:
            invalid += 1
            g = shapely.make_valid(g)
        g = scale(g, xfact=PX, yfact=PX, origin=(0, 0))
        if not g.is_valid:
            invalid += 1
            g = shapely.make_valid(g)
        if cl == "Islet":
            assert p["object_class"] in ["Islet", "EndocrineCluster"]
            cores[p["name"]] = g
        elif cl == "IsletExpanded":
            expanded[p["name"].removesuffix("_exp20um")] = g
        elif cl == "Tissue":
            tissue.append(g)
        else:
            parts[cl].append(g)
    assert set(cores) == set(expanded) == set(expected.object_name)
    tissue = shapely.union_all(tissue)
    names = list(cores)
    core_geoms = np.asarray([cores[n] for n in names], dtype=object)
    ctree = STRtree(core_geoms)
    arrays = {kind: np.asarray(v, dtype=object) for kind, v in parts.items()}
    trees = {kind: STRtree(v) for kind, v in arrays.items()}
    meta = expected.set_index("object_name")
    rows = []
    for name, core in cores.items():
        # Saved expanded geometry already encodes the app's tissue-clipped
        # 20 µm band. Explicitly subtract overlapping Islet/Cluster cores.
        # Singles are not subtracted, matching assemble_uniform.py's geoms.
        exp = expanded[name]
        ix = ctree.query(exp, predicate="intersects")
        local_cores = shapely.union_all(core_geoms[ix])
        peri = exp.difference(local_cores).intersection(tissue)
        m = meta.loc[name]
        assert np.isclose(core.area, m.core_area_um2, rtol=1e-7, atol=1e-5)
        for region, geom in [("core", core), ("peri", peri)]:
            denominator = float(geom.area)
            for kind in KINDS:
                near = trees[kind].query(geom, predicate="intersects")
                if len(near):
                    intersections = shapely.intersection(arrays[kind][near], geom)
                    areas = shapely.area(intersections)
                    positive = areas > 0
                    raw_sum = float(areas.sum())
                    union_area = float(shapely.union_all(intersections[positive]).area) if positive.any() else 0.
                    count = int(positive.sum())
                else:
                    raw_sum = union_area = 0.
                    count = 0
                assert union_area <= denominator + 1e-7 * max(denominator, 1)
                assert union_area >= 0
                rows.append({"islet_id": str(m.islet_id), "case_id": int(case),
                             "donor_status": str(m.donor_status), "object_name": name,
                             "object_class": str(m.object_class), "class": kind, "region": region,
                             "type": "islet_core" if region == "core" else "islet_band",
                             "core_area_um2": float(core.area), "area_um2": union_area,
                             "region_um2": denominator, "area_density": union_area / denominator if denominator > 0 else np.nan,
                             "area_percent": 100 * union_area / denominator if denominator > 0 else np.nan,
                             "raw_annotation_sum_area_um2": raw_sum,
                             "overlap_doublecount_area_um2": max(0., raw_sum - union_area),
                             "intersecting_annotation_count": count,
                             "core_complete": bool(m.core_complete),
                             "neighborhood_complete_20": bool(m.neighborhood_complete_20),
                             "paired_row_available": bool(m.paired_row_available),
                             "source_geometry": str(path)})
    result = pd.DataFrame(rows)
    audit = {"case_id": case, "objects": len(cores), "rows": len(result),
             "geometry_sha256": sha(path), "invalid_geometries_repaired_in_memory": invalid,
             "structure_annotation_counts": {k: len(v) for k, v in arrays.items()},
             "rows_with_overlap_area_gt_1e_6_um2": int((result.overlap_doublecount_area_um2 > 1e-6).sum()),
             "max_area_percent": float(result.area_percent.max())}
    return result, audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    census = pd.read_parquet(CENSUS)
    census = census[census.object_class.isin(["Islet", "EndocrineCluster"])].copy()
    assert len(census) == 4906 and census.islet_id.is_unique
    assert census.case_id.nunique() == 20 and not census.case_id.eq(78).any()
    pieces, audits = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(donor, int(case), frame): int(case) for case, frame in census.groupby("case_id")}
        for future in as_completed(futures):
            d, audit = future.result()
            pieces.append(d); audits.append(audit)
            print(f"Finished donor {audit['case_id']}: {audit['objects']} objects", flush=True)
    d = pd.concat(pieces, ignore_index=True).sort_values(["case_id", "islet_id", "region", "class"])
    assert len(d) == 4906 * 2 * 3
    assert not d.duplicated(["islet_id", "class", "region"]).any()
    assert d.area_percent.dropna().between(-1e-8, 100+1e-6).all()
    old = pd.read_parquet(OLD)
    old = old[old["class"].isin(KINDS) & old.type.isin(["islet_core", "islet_band"])].copy()
    old["islet_id"] = old["Case ID"].astype(int).astype(str) + "_" + old.islet_key
    old = old.rename(columns={"area_um2": "previous_summed_area_um2", "region_um2": "previous_region_um2"})
    d = d.merge(old[["islet_id", "class", "type", "previous_summed_area_um2", "previous_region_um2"]],
                on=["islet_id", "class", "type"], how="left", validate="one_to_one")
    assert d.previous_region_um2.notna().all()
    d["region_area_difference_um2"] = d.region_um2 - d.previous_region_um2
    d.to_csv(OUT / "structure_union_area_values.csv", index=False)
    pd.DataFrame(audits).sort_values("case_id").to_csv(OUT / "structure_union_donor_audit.csv", index=False)
    provenance = {"method": "union of intersected existing same-class annotations; clipped to exact saved core or saved20um expanded geometry minus all Islet/EndocrineCluster cores and intersected with Tissue; Singles are not subtracted, matching current app geometry",
                  "structure_classes": KINDS, "pixel_um": PX, "objects": 4906, "donors": 20, "rows": len(d),
                  "unpaired_objects_retained": int(census.paired_row_available.eq(False).sum()),
                  "max_absolute_region_area_difference_um2": float(d.region_area_difference_um2.abs().max()),
                  "rows_with_overlapping_annotation_area": int((d.overlap_doublecount_area_um2 > 1e-6).sum()),
                  "previous_rows_above_100_percent": int((100*d.previous_summed_area_um2/d.previous_region_um2 > 100+1e-6).sum()),
                  "corrected_rows_above_100_percent": int(d.area_percent.gt(100+1e-6).sum()),
                  "source_sha256": {str(p): sha(p) for p in [CENSUS, OLD, Path(__file__).resolve()]},
                  "donor_audits": sorted(audits, key=lambda x: x["case_id"])}
    (OUT / "structure_union_provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps({k: v for k, v in provenance.items() if k not in ["source_sha256", "donor_audits"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
