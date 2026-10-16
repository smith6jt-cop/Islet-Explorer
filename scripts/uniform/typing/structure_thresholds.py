#!/usr/bin/env python3
"""structure_thresholds.py -- structure objects for the five additions from thresholds fitted on the original 15.

Owner, 2026-10-16: "the pixel classifiers amount to a simple thresholder so for existing it is possible to find
optimal values to use for the new donors."

The original project's structure classifiers (`QuPath_project_IsletExplorer/classifiers/pixel_classifiers/`)
are linear models on smoothed intensity -- Nerve and Capillary a single-layer ANN_MLP (no hidden layer) on
Gaussian / Laplacian / weighted-SD / gradient features at sigma 0.5 and 1 px of B3TUBB+INS and CD31, Lymphatics a
50-tree RTrees on the same features of PDPN -- and `4_Islet_Cells.groovy` turned their output into annotations
inside Tissue (createAnnotationsFromPixelClassifier, minimum area 2 um2 for Nerve and 4 um2 for Capillary and
Lymphatics, SPLIT) and then dropped objects touching a 10 um ring at the image edge.

fit   On each ORIGINAL image, for each class, the threshold on the Gaussian-smoothed marker (sigma SIGMAS px at
      0.3774 um) that best reproduces that image's existing objects: maximum Dice of (smoothed >= t) against the
      objects' pixels, inside Tissue and outside the edge ring. The best sigma is the one with the highest median
      Dice. Written to `$OUT/structures/thresholds_fit.json`.
make  The thresholds agree across the original 15 in raw units, not after the UniFORM factors (the original
      classifiers were one fixed raw cut), so the median raw threshold of the 15 is the reference, moved into each
      ADDITION's intensity scale by the ratio of UniFORM factors (`$OUT/typing/shifts.json`):
      t_addition = t_median x median(factor over the 15) / factor(addition); and that image's objects are drawn as the original script drew them:
      smoothed >= t inside Tissue, 8-connected components split into separate objects, the class's minimum area,
      no object touching the 10 um edge ring. The additions' earlier objects (the same classifiers applied to
      their own raw intensities) are kept in `$OUT/structures/superseded_additions/` and compared.

    python scripts/uniform/typing/structure_thresholds.py --out $OUT fit
    python scripts/uniform/typing/structure_thresholds.py --out $OUT make
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import cv2
import numpy as np
import shapely
from scipy.ndimage import gaussian_filter
from shapely.geometry import Polygon, mapping, shape

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
from redsea_cells import Channels  # noqa: E402

CLASSES = {'Nerve': 'B3TUBB', 'Capillary': 'CD31', 'Lymphatic': 'PDPN'}
MIN_AREA_UM2 = {'Nerve': 2.0, 'Capillary': 4.0, 'Lymphatic': 4.0}
EDGE_MARGIN_UM = 10.0
SIGMAS = (0.5, 1.0, 1.5, 2.0, 3.0)
NBINS = 4000


def fill(mask: np.ndarray, geom: dict) -> None:
    """Scanline fill of a GeoJSON (Multi)Polygon into mask (1 inside, holes back to 0), at 1/16 px precision.
    skimage.draw.polygon tests every pixel of the bounding box against every vertex: hours for a Tissue outline."""
    polys = geom['coordinates'] if geom['type'] == 'MultiPolygon' else [geom['coordinates']]
    for poly in polys:
        rings = [np.round(np.asarray(r, dtype=np.float64) * 16).astype(np.int32) for r in poly if len(r) >= 3]
        if not rings:
            continue
        cv2.fillPoly(mask, [rings[0]], 1, lineType=cv2.LINE_8, shift=4)
        if len(rings) > 1:
            cv2.fillPoly(mask, rings[1:], 0, lineType=cv2.LINE_8, shift=4)


def tissue_mask(out: Path, tag: str, shape_):
    raw = json.loads((out / 'inputs' / 'tissue' / f'{tag}.geojson').read_text())
    feats = raw['features'] if isinstance(raw, dict) else raw
    m = np.zeros(shape_, dtype=np.uint8)
    for f in feats:
        if ((f.get('properties') or {}).get('classification') or {}).get('name') == 'Tissue':
            fill(m, f['geometry'])
    return m.astype(bool)


def edge_ring(shape_):
    H, W = shape_
    p = int(round(EDGE_MARGIN_UM / CT.PIXEL_UM))
    ring = np.ones(shape_, bool)
    ring[p:H - p, p:W - p] = False
    return ring


def object_mask(path: Path, shape_):
    m = np.zeros(shape_, dtype=np.uint8)
    for f in json.loads(path.read_text())['features']:
        if f.get('geometry'):
            fill(m, f['geometry'])
    return m.astype(bool)


def channel(out: Path, tag: str, marker: str):
    entry = CT.manifest(out)[tag]
    names = CT.channels(out, entry['entry_id'])
    img = Channels(entry['image_path'])
    try:
        return img(names.index(marker)), img.shape[1:]
    finally:
        img.close()


def best_threshold(values, inside):
    """Dice-optimal threshold of values against the boolean target, from two histograms."""
    lv = np.log1p(np.clip(values, 0, None))
    edges = np.linspace(0, float(lv.max()) + 1e-6, NBINS + 1)
    pos = np.histogram(lv[inside], bins=edges)[0][::-1].cumsum()[::-1]       # pixels >= edge, in target
    neg = np.histogram(lv[~inside], bins=edges)[0][::-1].cumsum()[::-1]
    total = inside.sum()
    dice = 2 * pos / np.maximum(2 * pos + neg + (total - pos), 1)
    k = int(np.argmax(dice))
    t = float(np.expm1(edges[k]))
    prec = pos[k] / max(pos[k] + neg[k], 1)
    rec = pos[k] / max(total, 1)
    return t, float(dice[k]), float(prec), float(rec)


def fit_one(out: Path, tag: str) -> dict:
    res = {}
    for cls, marker in CLASSES.items():
        ch, shp = channel(out, tag, marker)
        region = tissue_mask(out, tag, shp) & ~edge_ring(shp)
        target = object_mask(out / 'structures' / cls / f'{cls}__{tag}.geojson', shp)[region]
        res[cls] = {'marker': marker, 'object_pixels': int(target.sum())}
        for s in SIGMAS:
            sm = gaussian_filter(ch, s)[region]
            t, d, p, r = best_threshold(sm, target)
            res[cls][str(s)] = {'threshold': t, 'dice': d, 'precision': p, 'recall': r}
        del ch
    print(tag, json.dumps({c: {s: round(v[s]['dice'], 3) for s in map(str, SIGMAS)} for c, v in res.items()}),
          flush=True)
    return res


def polygons(binary: np.ndarray, min_area_px: float) -> list:
    """8-connected components of a binary image as polygons with holes, pixel coordinates."""
    contours, hier = cv2.findContours(binary.astype(np.uint8), cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    out = []
    if hier is None:
        return out
    hier = hier[0]
    for i, c in enumerate(contours):
        if hier[i][3] != -1:            # a hole; handled with its parent
            continue
        shell = c[:, 0, :].astype(float)
        if len(shell) < 3:
            continue
        holes = []
        j = hier[i][2]
        while j != -1:
            h = contours[j][:, 0, :].astype(float)
            if len(h) >= 3:
                holes.append(h)
            j = hier[j][0]
        g = Polygon(shell, holes).buffer(0)
        if g.is_empty or g.area < min_area_px:
            continue
        out.append(g)
    return out


def make_one(out: Path, tag: str, thresholds: dict, sigma: float) -> dict:
    rec = {}
    for cls, marker in CLASSES.items():
        ch, shp = channel(out, tag, marker)
        tissue = tissue_mask(out, tag, shp)
        binary = (gaussian_filter(ch, sigma) >= thresholds[cls]) & tissue
        del ch
        geoms = polygons(binary, MIN_AREA_UM2[cls] / CT.PIXEL_UM ** 2)
        H, W = shp
        p = EDGE_MARGIN_UM / CT.PIXEL_UM
        inner = shapely.box(p, p, W - p, H - p)
        geoms = [g for g in geoms if inner.contains(g)]
        feats = [{'type': 'Feature', 'geometry': mapping(g),
                  'properties': {'objectType': 'annotation', 'classification': {'name': cls}}} for g in geoms]
        dest = out / 'structures' / cls / f'{cls}__{tag}.geojson'
        old = out / 'structures' / 'superseded_additions' / cls / f'{cls}__{tag}.geojson'
        if not old.exists():
            old.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(dest, old)
        draw_mask = np.zeros(shp, np.uint8)
        for f in feats:
            fill(draw_mask, f['geometry'])
        new_mask = draw_mask.astype(bool)
        old_mask = object_mask(old, shp)
        inter = (new_mask & old_mask).sum()
        rec[cls] = {'threshold_raw': thresholds[cls], 'objects': len(feats),
                    'objects_before': len(json.loads(old.read_text())['features']),
                    'area_um2': float(new_mask.sum() * CT.PIXEL_UM ** 2),
                    'area_um2_before': float(old_mask.sum() * CT.PIXEL_UM ** 2),
                    'dice_vs_before': float(2 * inter / max(new_mask.sum() + old_mask.sum(), 1))}
        dest.write_text(json.dumps({'type': 'FeatureCollection', 'features': feats}))
    print(tag, json.dumps(rec), flush=True)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('mode', choices=['fit', 'make'])
    ap.add_argument('--workers', type=int, default=5)
    a = ap.parse_args()
    rows = CT.cohort_rows()
    originals = [r['image_name'] for r in rows if r['cohort'] == 'original']
    additions = [r['image_name'] for r in rows if r['cohort'] == 'addition']
    root = a.out / 'structures'
    if a.mode == 'fit':
        with ProcessPoolExecutor(a.workers) as pool:
            fits = dict(zip(originals, pool.map(fit_one, [a.out] * len(originals), originals)))
        sig = {str(s): float(np.median([fits[t][c][str(s)]['dice'] for t in originals for c in CLASSES]))
               for s in SIGMAS}
        (root / 'thresholds_fit.json').write_text(json.dumps({'sigmas_median_dice': sig, 'images': fits}, indent=1))
        print('median Dice by sigma:', sig)
        return 0
    fit = json.loads((root / 'thresholds_fit.json').read_text())
    sigma = max(fit['sigmas_median_dice'], key=fit['sigmas_median_dice'].get)
    shifts = json.loads((CT.typing_dir(a.out) / 'shifts.json').read_text())['shifts']
    common, per_image = {}, {}
    def spread(v):
        v = np.asarray(list(v))
        return float(np.percentile(v, 90) / np.percentile(v, 10))
    for cls, marker in CLASSES.items():
        f = shifts[marker]['scale']
        raw = {t: fit['images'][t][cls][sigma]['threshold'] for t in originals}
        aligned = {t: raw[t] * f[t] for t in originals}
        # Measured 2026-10-16: the fitted thresholds agree across the original 15 in RAW units (p90/p10 1.4 Nerve,
        # 1.4 Capillary, 2.0 Lymphatic) and not after the UniFORM factors (4.6, 2.8, 14): the original classifiers
        # were one fixed raw cut on images acquired alike. So the originals' raw scale is the reference, and an
        # addition's threshold is that raw median moved by the addition's intensity scale relative to the
        # originals' (the ratio of UniFORM factors).
        f_ref = float(np.median([f[t] for t in originals]))
        common[cls] = float(np.median(list(raw.values())))
        per_image[cls] = {'raw_by_original': raw, 'raw_spread_p90_over_p10': spread(raw.values()),
                          'aligned_by_original': aligned, 'aligned_spread_p90_over_p10': spread(aligned.values()),
                          'originals_median_factor': f_ref,
                          'raw_by_addition': {t: common[cls] * f_ref / f[t] for t in additions}}
    made = {}
    with ProcessPoolExecutor(a.workers) as pool:
        futs = {t: pool.submit(make_one, a.out, t, {c: per_image[c]['raw_by_addition'][t] for c in CLASSES},
                               float(sigma)) for t in additions}
        for t, fu in futs.items():
            made[t] = fu.result()
    (root / 'thresholded.json').write_text(json.dumps({
        'sigma_px': float(sigma), 'min_area_um2': MIN_AREA_UM2, 'edge_margin_um': EDGE_MARGIN_UM,
        'reference_raw_threshold': common, 'thresholds': per_image, 'images': made}, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
