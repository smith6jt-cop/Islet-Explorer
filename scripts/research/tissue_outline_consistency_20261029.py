#!/usr/bin/env python3
"""Are the 20 images' Tissue outlines drawn the same way? (owner 2026-10-29: "the tissue outlines need to be consistent
across the images")

Each image's outline in use (`<run>/inputs/tissue/<image>.geojson`) is compared with one rule applied to every image,
redraw_tissue.py's (pyramid level 3, nuclei from DAPI or stroma from ColIV/VIM/SMA, holes under 1,000 um2 filled,
pieces under 80 um2 removed), once keeping the larger holes (as the original 15 outlines do) and once filling every hole
(as the five additions' redrawn outlines do). Per image: Dice with each version; the holes of the outline in use (count,
area as a proportion of the outline with its holes filled) and how much of the hole area holds nuclei or stroma by the
same rule (tissue the outline leaves out, where no cell was detected); the rule's tissue outside the outline in use.

Writes output/tissue_outline_consistency_20261029/consistency.csv and <image>.png (DAPI blue, stroma green; outline in
use white; its holes with nuclei red, with stroma only yellow).

    python scripts/research/tissue_outline_consistency_20261029.py
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage
from skimage import filters, morphology

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform'))
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402
import redraw_tissue as RT  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'tissue_outline_consistency_20261029'


def rule_parts(dapi: np.ndarray, stroma: np.ndarray):
    d = np.log1p(filters.gaussian(dapi.astype(float), 2, preserve_range=True))
    s = filters.gaussian(stroma, 2)
    nuclei, strom = RT.above(d, d[d > 0], 0.6), RT.above(s, s, 1.2)
    m = morphology.remove_small_holes(nuclei | strom, int(1000 / RT.PX2))
    m = morphology.remove_small_objects(m, int(RT.MIN_PIECE_UM2 / RT.PX2))
    return nuclei, strom, m, ndimage.binary_fill_holes(m)


def outline_mask(tag: str, shape: tuple[int, int]) -> np.ndarray:
    raw = json.loads((RUN / 'inputs' / 'tissue' / f'{tag}.geojson').read_text())
    feats = raw['features'] if isinstance(raw, dict) else raw
    m = np.zeros(shape, np.uint8)
    for f in feats:
        if ((f.get('properties') or {}).get('classification') or {}).get('name') != 'Tissue':
            continue
        g = f['geometry']
        for poly in (g['coordinates'] if g['type'] == 'MultiPolygon' else [g['coordinates']]):
            rings = [np.round(np.asarray(r, float) / RT.DOWNSAMPLE * 16).astype(np.int32) for r in poly if len(r) >= 3]
            if rings:
                cv2.fillPoly(m, [rings[0]], 1, lineType=cv2.LINE_8, shift=4)
                if len(rings) > 1:
                    cv2.fillPoly(m, rings[1:], 0, lineType=cv2.LINE_8, shift=4)
    return m.astype(bool)


def dice(a: np.ndarray, b: np.ndarray) -> float:
    return float(2 * (a & b).sum() / max(a.sum() + b.sum(), 1))


def one(tag: str) -> dict:
    dapi, stroma, _ = RT.channels(RT.IMAGES / f'{tag}.ome.tiff')
    nuclei, strom, keep, fill = rule_parts(dapi, stroma)
    cur = outline_mask(tag, dapi.shape)
    env = ndimage.binary_fill_holes(cur)
    holes = env & ~cur
    lab, n = ndimage.label(holes)
    sizes = np.bincount(lab.ravel())[1:] * RT.PX2
    big = np.isin(lab, np.flatnonzero(sizes >= 1000) + 1)
    rec = {'image': tag, 'outline_mm2': cur.sum() * RT.PX2 / 1e6,
           'dice_rule_holes_kept': dice(cur, keep), 'dice_rule_holes_filled': dice(cur, fill),
           'holes_ge_1000um2': int((sizes >= 1000).sum()), 'holes_lt_1000um2': int((sizes < 1000).sum()),
           'hole_area_pct': 100 * holes.sum() / max(env.sum(), 1),
           'hole_area_with_nuclei_pct': 100 * (holes & nuclei).sum() / max(env.sum(), 1),
           'hole_area_with_stroma_only_pct': 100 * (holes & strom & ~nuclei).sum() / max(env.sum(), 1),
           'rule_tissue_outside_outline_pct': 100 * (keep & ~env).sum() / max(env.sum(), 1),
           'outline_outside_rule_pct': 100 * (cur & ~fill).sum() / max(env.sum(), 1)}
    OUT.mkdir(parents=True, exist_ok=True)
    def u8(a):
        lo, hi = np.percentile(a, [1, 99.5])
        return (np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)
    rgb = np.dstack([np.zeros(dapi.shape, np.uint8), u8(stroma), u8(np.log1p(dapi.astype(float)))])
    over = rgb.copy()
    over[holes & nuclei] = (255, 0, 0)
    over[holes & strom & ~nuclei] = (255, 255, 0)
    edge = cv2.morphologyEx(cur.astype(np.uint8), cv2.MORPH_GRADIENT, np.ones((3, 3), np.uint8)).astype(bool)
    over[edge] = (255, 255, 255)
    cv2.imwrite(str(OUT / f'{tag}.png'), np.hstack([rgb, over])[:, :, ::-1])
    print(tag, {k: round(v, 3) if isinstance(v, float) else v for k, v in rec.items() if k != 'image'}, flush=True)
    return rec


def main() -> int:
    tags = [r['image_name'] for r in CT.cohort_rows()]
    with ProcessPoolExecutor(5) as ex:
        rows = list(ex.map(one, tags))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / 'consistency.csv', index=False)
    pd.set_option('display.width', 250)
    print(d.round(3).to_string(index=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
