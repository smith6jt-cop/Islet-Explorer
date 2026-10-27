#!/usr/bin/env python3
"""Redraw the Tissue outline of the five added donors' images (owner, 2026-10-27).

The additions' Tissue outlines came from the `Tissue` pixel classifier of the six-donor review project and were reused
unchanged; ND_6232's covered 38 % of the tissue visible in DAPI, and ND_6539 / Aab_6558 held thousands of cell-sized
holes (CLAUDE.md, 2026-10-26 standing rule). The redraw, approved by the owner "if all holes are filled":

  - pyramid level 3 (3.0192 um/px, the original Tissue classifier's resolution);
  - nuclei: DAPI, Gaussian sigma 2 px, log; Otsu on the nonzero pixels, threshold moved to 0.6 of the way from the
    5th percentile to Otsu's;
  - stroma: ColIV, VIM and SMA, each log-scaled to its own 1st-99th percentiles in the image, averaged, Gaussian
    sigma 2 px; Otsu threshold at 1.2 of the way from the 5th percentile;
  - Tissue = nuclei OR stroma; pieces under 80 um2 removed (the original 15 outlines' smallest piece is 82 um2); every
    hole filled (the original 15 keep holes of 1,000 um2 and more; the owner's rule for the redraw fills all);
  - clipped to the image less a 5 um border, as the original Tissue was (`edgeMarginUm` and `clipObjectsToInner` in
    the original project's 1_Panc_Objects.groovy, "to fully cover the bright rim").

Without the hole filling the same rule reproduces the original 15 outlines at Dice 0.958-0.991.

Writes <out-dir>/<image>.geojson: one Tissue annotation (MultiPolygon, full-resolution pixel coordinates) that
run_uniform.py's tissue_export step copies for images whose cohort.csv `tissue_outline` names it.

    python scripts/uniform/redraw_tissue.py --images Aab_6090 ND_6232 T1D_6520 ND_6539 Aab_6558 \\
        --out-dir scripts/uniform/tissue_outlines
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import cv2
import numpy as np
import tifffile
from scipy import ndimage
from shapely.geometry import MultiPolygon, Polygon, box, mapping
from shapely.ops import unary_union
from shapely.validation import make_valid
from skimage import filters, morphology

REPO = Path(__file__).resolve().parents[2]
IMAGES = REPO / 'data/development/uniform-20donor-20261022-v10/images'
LEVEL, DOWNSAMPLE = 3, 8
PIXEL_UM = 0.3774
PX2 = (PIXEL_UM * DOWNSAMPLE) ** 2
STROMA = ('ColIV', 'VIM', 'SMA')
MIN_PIECE_UM2 = 80
EDGE_MARGIN_UM = 5.0


def channels(path: Path) -> tuple[np.ndarray, np.ndarray, tuple[int, int]]:
    with tifffile.TiffFile(path) as t:
        height, width = t.series[0].levels[0].shape[-2:]
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        level = t.series[0].levels[LEVEL]
        if level.shape[-2:] != tuple(-(-s // DOWNSAMPLE) for s in t.series[0].levels[0].shape[-2:]) and \
                level.shape[-2:] != tuple(s // DOWNSAMPLE for s in t.series[0].levels[0].shape[-2:]):
            raise SystemExit(f'{path}: pyramid level {LEVEL} is not a {DOWNSAMPLE}x downsample')
        dapi = level.pages[names.index('DAPI')].asarray()
        stroma = []
        for n in STROMA:
            v = np.log1p(level.pages[names.index(n)].asarray().astype(np.float32))
            lo, hi = np.percentile(v, [1, 99])
            stroma.append(np.clip((v - lo) / (hi - lo + 1e-6), 0, 1))
    return dapi, np.mean(stroma, axis=0), (width, height)


def above(v: np.ndarray, fg: np.ndarray, k: float) -> np.ndarray:
    t, bg = filters.threshold_otsu(fg), np.percentile(fg, 5)
    return v > bg + k * (t - bg)


def tissue_mask(dapi: np.ndarray, stroma: np.ndarray) -> np.ndarray:
    d = np.log1p(filters.gaussian(dapi.astype(float), 2, preserve_range=True))
    s = filters.gaussian(stroma, 2)
    m = above(d, d[d > 0], 0.6) | above(s, s, 1.2)
    # the order the owner reviewed: holes under 1,000 um2 filled, small pieces removed, then every hole filled
    m = morphology.remove_small_holes(m, int(1000 / PX2))
    m = morphology.remove_small_objects(m, int(MIN_PIECE_UM2 / PX2))
    return ndimage.binary_fill_holes(m)


def to_geometry(mask: np.ndarray, size: tuple[int, int]) -> MultiPolygon:
    """Outer contours traced at 4x (OpenCV traces pixel centres: a quarter-pixel error), corner contacts repaired."""
    scale = DOWNSAMPLE / 4
    cs, hier = cv2.findContours(np.kron(mask.astype(np.uint8), np.ones((4, 4), np.uint8)), cv2.RETR_CCOMP,
                                cv2.CHAIN_APPROX_SIMPLE)
    parts = []
    for i, c in enumerate(cs):
        if hier[0][i][3] != -1 or len(c) < 3:
            continue
        holes, j = [], hier[0][i][2]
        while j != -1:
            if len(cs[j]) >= 3:
                holes.append((cs[j][:, 0, :] + 0.5) * scale)
            j = hier[0][j][0]
        parts.append(make_valid(Polygon((c[:, 0, :] + 0.5) * scale, holes)))
    m = EDGE_MARGIN_UM / PIXEL_UM
    geom = unary_union(parts).intersection(box(m, m, size[0] - m, size[1] - m))
    geom = MultiPolygon([p for g in getattr(geom, 'geoms', [geom]) for p in getattr(g, 'geoms', [g])
                         if p.geom_type == 'Polygon'])
    if not geom.is_valid:
        raise SystemExit('outline geometry is not valid')
    return geom


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--out-dir', type=Path, required=True)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    for tag in a.images:
        dapi, stroma, size = channels(IMAGES / f'{tag}.ome.tiff')
        mask = tissue_mask(dapi, stroma)
        geom = to_geometry(mask, size)
        holes = sum(len(p.interiors) for p in geom.geoms)
        feature = {'type': 'Feature', 'geometry': mapping(geom),
                   # exactly the original exports' properties: an unnamed annotation, so QuPath's cell `Parent` reads
                   # "Annotation (Tissue)", the value the pipeline keys tissue cells on (islet_cell_test.TISSUE_PARENT)
                   'properties': {'objectType': 'annotation', 'classification': {'name': 'Tissue', 'color': [102, 153, 153]},
                                  'isLocked': True}}
        out = a.out_dir / f'{tag}.geojson'
        out.write_text(json.dumps({'type': 'FeatureCollection', 'features': [feature]}))
        print(f'{tag}: {mask.sum() * PX2 / 1e6:.2f} mm2 (mask), {geom.area * PIXEL_UM ** 2 / 1e6:.2f} mm2 (outline), '
              f'{len(geom.geoms)} pieces, {holes} holes -> {out}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
