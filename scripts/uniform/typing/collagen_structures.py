#!/usr/bin/env python3
"""Collagen IV structure objects, drawn the way structure_thresholds.py draws Nerve and Capillary.

Owner 2026-10-29: "ColIV needs to be converted to objects like vessels, lymphatics, and nerves." There is no original
classifier for ColIV to fit a threshold against, so each image sets its own, by the typing's known-negative method
(negative_p99): the 99th percentile of the smoothed channel over the pixels of isolated Exocrine cells (typing/identity),
at least ISOLATION_UM from every Capillary, Lymphatic and Nerve object and every endocrine object. Every acinus is
outlined by a dim collagen IV basement membrane; an object is collagen IV brighter than that mesh. Agent determination
2026-10-29, compared on ND_0112, ND_6232 and T1D_6533 against the background peak plus 1, 2 or 3 widths (3 widths drew
0.2-2 % of tissue and missed most vessels; 1 width outlined the acinar mesh; 2 widths reached 44-93 % of capillary
pixels), a Gaussian-mixture or Otsu split (42-66 % of tissue, one tissue-wide network) and the triangle method
(2-95 %). Drawing as structure_thresholds.make: the despeckled channel (redsea_cells.despeckle; owner 2026-10-28)
Gaussian-smoothed at SIGMA px (the 15's best Dice for Nerve and Capillary), >= threshold inside Tissue, 8-connected
components as separate objects, at least MIN_AREA_UM2 (Capillary's), none touching the 10 um edge ring.

Writes structures/Collagen IV/Collagen IV__<image>.geojson, structures/collagen_iv.json (threshold and checks per
image) and output/structure_check_20261029/collagen_iv/ renders (whole image and zoomed tiles: ColIV with DAPI and
the object outlines) for the image check every gating geometry needs (CLAUDE.md, 2026-10-26).

    python scripts/uniform/typing/collagen_structures.py --run data/development/uniform-20donor-20261027 --workers 4
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import shapely
from scipy.ndimage import gaussian_filter
from shapely.geometry import mapping

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
import redsea_cells as RS  # noqa: E402
import registry as RG  # noqa: E402
import structure_thresholds as ST  # noqa: E402

CLASS = 'Collagen IV'
MARKER = 'ColIV'
SIGMA = 2.0
MIN_AREA_UM2 = 4.0
REFERENCE_QUANTILE = 99.0           # negative_p99
REFERENCE_STRUCTURES = ('Capillary', 'Lymphatic', 'Nerve')
REPO = Path(__file__).resolve().parents[3]
RENDERS = REPO / 'output' / 'structure_check_20261029' / 'collagen_iv'


def load(run: Path, tag: str, marker: str) -> np.ndarray:
    entry = CT.manifest(run)[tag]
    names = CT.channels(run, entry['entry_id'])
    img = RS.Channels(entry['image_path'])
    try:
        return RS.despeckle(img.native(names.index(marker))).astype(np.float32)
    finally:
        img.close()


def stretch(a: np.ndarray, region: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(a[region], [1, 99.5]) if region.any() else (0, 1)
    return (np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)


def render(tag: str, col: np.ndarray, dapi: np.ndarray, tissue: np.ndarray, binary: np.ndarray, capillary: np.ndarray):
    """Whole image (downsampled) and four zoomed tiles: ColIV green, DAPI blue, Capillary objects' outline yellow, the
    Collagen IV objects' outline magenta on top."""
    RENDERS.mkdir(parents=True, exist_ok=True)
    g, b = stretch(col, tissue), stretch(dapi, tissue)
    rgb = np.dstack([np.zeros_like(g), g, b])
    edge = cv2.morphologyEx(binary.astype(np.uint8), cv2.MORPH_GRADIENT, np.ones((3, 3), np.uint8)).astype(bool)
    cap = cv2.morphologyEx(capillary.astype(np.uint8), cv2.MORPH_GRADIENT, np.ones((3, 3), np.uint8)).astype(bool)
    over = rgb.copy()
    over[cap] = (255, 255, 0)
    over[edge] = (255, 0, 255)
    H, W = col.shape
    f = 2000 / max(H, W)
    small = cv2.resize(np.hstack([rgb, over]), (int(2 * W * f), int(H * f)), interpolation=cv2.INTER_AREA)
    cv2.imwrite(str(RENDERS / f'{tag}_whole.png'), small[:, :, ::-1])
    ys, xs = np.nonzero(binary[::50, ::50] & tissue[::50, ::50])
    rng = np.random.default_rng(0)
    tiles = []
    for i in rng.choice(len(ys), size=min(4, len(ys)), replace=False) if len(ys) else []:
        y, x = min(ys[i] * 50, H - 800), min(xs[i] * 50, W - 800)
        tiles.append(np.hstack([rgb[y:y + 800, x:x + 800], over[y:y + 800, x:x + 800]]))
    if tiles:
        cv2.imwrite(str(RENDERS / f'{tag}_tiles.png'), np.vstack(tiles)[:, :, ::-1])


def reference(run: Path, tag: str, shp, tissue: np.ndarray) -> np.ndarray:
    """Pixels of the image's isolated Exocrine cells: at least ISOLATION_UM from every structure and endocrine object."""
    iso = int(round(RG.load().const('ISOLATION_UM') / CT.PIXEL_UM))
    feats = json.loads((run / 'cell_geojson' / f'{tag}.geojson').read_text())
    feats = feats['features'] if isinstance(feats, dict) else feats
    mask, _, oids, _ = RS.rasterize(feats, shp)
    del feats
    ident = pd.read_parquet(run / 'typing' / 'identity' / f'{tag}.parquet', columns=['object_id', 'label'])
    exocrine = set(ident.loc[ident['label'] == 'Exocrine', 'object_id'])
    keep = np.zeros(len(oids) + 1, bool)
    keep[[i + 1 for i, o in enumerate(oids) if o in exocrine]] = True
    near = np.zeros(shp, np.uint8)
    for cls in REFERENCE_STRUCTURES:
        near |= ST.object_mask(run / 'structures' / cls / f'{cls}__{tag}.geojson', shp).astype(np.uint8)
    raw = json.loads((run / 'islets' / f'{tag}.geojson').read_text())
    for f in raw['features'] if isinstance(raw, dict) else raw:
        if f.get('geometry'):
            ST.fill(near, f['geometry'])
    near = cv2.dilate(near, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * iso + 1, 2 * iso + 1))).astype(bool)
    return keep[mask] & ~near & tissue


def one(run: Path, tag: str) -> dict:
    col = load(run, tag, MARKER)
    shp = col.shape
    tissue = ST.tissue_mask(run, tag, shp)
    sm = gaussian_filter(col, SIGMA)
    ref = reference(run, tag, shp, tissue)
    if ref.sum() < 10000:
        raise ValueError(f'{tag}: {int(ref.sum())} reference pixels')
    threshold = float(np.percentile(sm[ref], REFERENCE_QUANTILE))
    binary = (sm >= threshold) & tissue
    del sm
    geoms = ST.polygons(binary, MIN_AREA_UM2 / CT.PIXEL_UM ** 2)
    H, W = shp
    p = ST.EDGE_MARGIN_UM / CT.PIXEL_UM
    inner = shapely.box(p, p, W - p, H - p)
    geoms = [g for g in geoms if inner.contains(g)]
    feats = [{'type': 'Feature', 'geometry': mapping(g),
              'properties': {'objectType': 'annotation', 'classification': {'name': CLASS}}} for g in geoms]
    dest = run / 'structures' / CLASS / f'{CLASS}__{tag}.geojson'
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps({'type': 'FeatureCollection', 'features': feats}))
    drawn = np.zeros(shp, np.uint8)
    for f in feats:
        ST.fill(drawn, f['geometry'])
    drawn = drawn.astype(bool)
    capillary = ST.object_mask(run / 'structures' / 'Capillary' / f'Capillary__{tag}.geojson', shp)
    near = cv2.dilate(drawn.astype(np.uint8), np.ones((7, 7), np.uint8)).astype(bool)   # within ~1 um
    areas = np.array([g.area for g in geoms]) * CT.PIXEL_UM ** 2
    rec = {'threshold_raw': threshold, 'reference_pixels': int(ref.sum()), 'objects': len(feats), 'tissue_pct': float(100 * (drawn & tissue).sum() / tissue.sum()),
           'object_area_um2_median': float(np.median(areas)) if len(areas) else np.nan,
           'object_area_um2_p99': float(np.percentile(areas, 99)) if len(areas) else np.nan,
           'largest_object_pct_of_objects': float(100 * areas.max() / areas.sum()) if len(areas) else np.nan,
           'capillary_pixels_within_1um_pct': float(100 * (capillary & near).sum() / max(capillary.sum(), 1)),
           'objects_per_mm2': float(len(feats) / (tissue.sum() * CT.PIXEL_UM ** 2 / 1e6))}
    render(tag, col, load(run, tag, 'DAPI'), tissue, drawn, capillary)
    print(tag, json.dumps({k2: round(x, 3) if isinstance(x, float) else x for k2, x in rec.items()}), flush=True)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--images', nargs='*', default=[])
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    tags = a.images or [r['image_name'] for r in CT.cohort_rows()]
    with ProcessPoolExecutor(a.workers) as ex:
        res = dict(zip(tags, ex.map(one, repeat(a.run), tags)))
    path = a.run / 'structures' / 'collagen_iv.json'
    images = json.loads(path.read_text())['images'] if path.exists() else {}
    images.update(res)
    path.write_text(json.dumps({'method': __doc__.split('\n\n')[1], 'sigma_px': SIGMA, 'min_area_um2': MIN_AREA_UM2,
                                'images': images}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
