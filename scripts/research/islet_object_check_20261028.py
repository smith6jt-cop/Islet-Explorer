#!/usr/bin/env python3
"""Check the endocrine objects that changed between the owner-reviewed run and the current run, against the image.

Owner standing rule (2026-10-26): geometry that gates data is checked against the image before it is used. The owner
reviewed the endocrine objects of run uniform-20donor-20261015 on 2026-10-16 ("The islet annotations look good"); the
2026-10-27 rerun (run uniform-20donor-20261027, after the five additions' Tissue outlines were redrawn) rewrote them.
For each image named, this compares the two object sets:

- metrics (both versions, same definitions): objects per class, total area, Dice between the two unions; hormone-high
  cells = REDSEA whole-cell max over INS/GCG/SST of value / that hormone's p99 among the image's cells outside every
  object, at >= 2; the proportion of hormone-high cells inside an object; objects holding no hormone-high cell;
- renders (output/islet_object_check_20261028/): <image>_whole.png, the whole image at 1.51 um/px, and
  <image>_tiles.png, full-resolution 250 um tiles at the six largest differences (object area in one version only).
  DAPI dim gray, insulin green, glucagon magenta, somatostatin blue, each windowed from the image's tissue median
  (black) to its p99.9 (white), gamma 1; outlines: current run yellow, reviewed run cyan.

    python scripts/research/islet_object_check_20261028.py --images Aab_6558 ND_6232 --workers 6
"""
from __future__ import annotations

import argparse
import json
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import tifffile
from PIL import Image, ImageDraw
from shapely import STRtree
from shapely.geometry import Point, shape, box
from shapely.ops import unary_union

REPO = Path(__file__).resolve().parents[2]
NEW = REPO / 'data/development/uniform-20donor-20261027'
OLD = REPO / 'data/development/uniform-20donor-20261015'
IMAGES = REPO / 'data/development/uniform-20donor-20261022-v10/images'
OUTDIR = REPO / 'output/islet_object_check_20261028'
PIXEL_UM = 0.3774
CLASSES = ('Islet', 'EndocrineCluster', 'EndocrineSingle')
HORMONES = ('INS', 'GCG', 'SST')
COLORS = {'DAPI': (0.30, 0.30, 0.30), 'INS': (0.0, 1.0, 0.0), 'GCG': (1.0, 0.0, 0.85), 'SST': (0.25, 0.55, 1.0)}
NEW_COLOR, OLD_COLOR = (255, 230, 0), (0, 230, 255)
TILE_UM = 250.0


def objects(run: Path, tag: str) -> list[tuple[str, object]]:
    raw = json.loads((run / 'islets' / f'{tag}.geojson').read_text())
    out = []
    for f in (raw['features'] if isinstance(raw, dict) else raw):
        c = f.get('properties', {}).get('classification', {})
        name = c.get('name') if isinstance(c, dict) else c
        if name in CLASSES:
            out.append((name, shape(f['geometry'])))
    return out


def cells(run: Path, tag: str) -> pd.DataFrame:
    xy = pq.read_table(run / 'cells' / f'{tag}.parquet',
                       columns=['Object ID', 'Centroid X µm', 'Centroid Y µm']).to_pandas()
    xy.columns = ['object_id', 'x', 'y']
    rs = pq.read_table(run / 'typing' / 'redsea' / f'{tag}.parquet',
                       columns=['object_id', 'qc_analysis_eligible'] + [f'Cell__{h}' for h in HORMONES]).to_pandas()
    d = xy.merge(rs[rs['qc_analysis_eligible']], on='object_id')
    d['px'], d['py'] = d['x'] / PIXEL_UM, d['y'] / PIXEL_UM
    return d


def membership(d: pd.DataFrame, objs: list) -> np.ndarray:
    """Index of the object holding each cell's centroid, -1 outside every object."""
    tree = STRtree([g for _, g in objs])
    pts = [Point(x, y) for x, y in zip(d['px'], d['py'])]
    idx = np.full(len(d), -1)
    hit_pt, hit_obj = tree.query(pts, predicate='within')
    idx[hit_pt] = hit_obj
    return idx


def metrics(tag: str, run: Path, objs: list) -> dict:
    d = cells(run, tag)
    idx = membership(d, objs)
    outside = idx < 0
    stat = np.zeros(len(d))
    for h in HORMONES:
        v = d[f'Cell__{h}'].to_numpy(float)
        p99 = np.nanpercentile(v[outside], 99)
        stat = np.maximum(stat, np.nan_to_num(v / p99 if p99 > 0 else v * 0))
    high = stat >= 2
    held = np.bincount(idx[high & ~outside], minlength=len(objs))
    r = {'cells': len(d), 'hormone_high_cells': int(high.sum()),
         'hormone_high_inside_objects': float(np.mean(~outside[high])) if high.any() else np.nan}
    for cls in CLASSES:
        sel = [i for i, (n, _) in enumerate(objs) if n == cls]
        r[f'{cls}_n'] = len(sel)
        r[f'{cls}_without_hormone_high'] = int(sum(held[i] == 0 for i in sel))
    r['area_um2'] = float(sum(g.area for _, g in objs) * PIXEL_UM ** 2)
    return r


def window(a: np.ndarray, tissue: np.ndarray) -> np.ndarray:
    lo = float(np.median(a[tissue])) if tissue.any() else float(np.median(a))
    hi = float(np.percentile(a, 99.9))
    return np.clip((a.astype(np.float32) - lo) / max(hi - lo, 1.0), 0, 1)


def composite(planes: dict, tissue: np.ndarray) -> np.ndarray:
    rgb = np.zeros(planes['DAPI'].shape + (3,), np.float32)
    for ch, col in COLORS.items():
        w = window(planes[ch], tissue)
        for k in range(3):
            rgb[..., k] += w * col[k]
    return (np.clip(rgb, 0, 1) * 255).astype(np.uint8)


def draw(img: Image.Image, objs: list, scale: float, offset=(0, 0), color=(255, 255, 0), width=2):
    dr = ImageDraw.Draw(img)
    for _, g in objs:
        for p in getattr(g, 'geoms', [g]):
            xy = [((x - offset[0]) * scale, (y - offset[1]) * scale) for x, y in p.exterior.coords]
            if len(xy) > 2:
                dr.line(xy + [xy[0]], fill=color, width=width)


def render(tag: str, new: list, old: list) -> None:
    with tifffile.TiffFile(IMAGES / f'{tag}.ome.tiff') as t:
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        lv2, lv0 = t.series[0].levels[2], t.series[0].levels[0]
        H0, W0 = lv0.shape[-2:]
        H2, W2 = lv2.shape[-2:]
        planes2 = {ch: lv2.pages[names.index(ch)].asarray() for ch in COLORS}
        tissue2 = planes2['DAPI'] > np.percentile(planes2['DAPI'], 50)
        img = Image.fromarray(composite(planes2, tissue2))
        s = W2 / W0
        draw(img, old, s, color=OLD_COLOR, width=2)
        draw(img, new, s, color=NEW_COLOR, width=2)
        img.save(OUTDIR / f'{tag}_whole.png')
        # the six largest differences: area of one version not covered by the other
        un, uo = unary_union([g for _, g in new]), unary_union([g for _, g in old])
        diff = [(g.area, 'current only', g) for g in getattr(un.difference(uo), 'geoms', [un.difference(uo)])] + \
               [(g.area, 'reviewed only', g) for g in getattr(uo.difference(un), 'geoms', [uo.difference(un)])]
        diff = sorted([x for x in diff if not x[2].is_empty], key=lambda x: -x[0])[:6]
        half = TILE_UM / PIXEL_UM / 2
        tiles = []
        full = {ch: lv0.pages[names.index(ch)] for ch in COLORS}
        for area, kind, g in diff:
            cx, cy = g.centroid.x, g.centroid.y
            x0, y0 = int(max(0, cx - half)), int(max(0, cy - half))
            x1, y1 = int(min(W0, cx + half)), int(min(H0, cy + half))
            planes = {ch: full[ch].asarray()[y0:y1, x0:x1] for ch in COLORS}
            tissue = planes['DAPI'] > np.percentile(planes['DAPI'], 50)
            tile = Image.fromarray(composite(planes, tissue))
            win = box(x0, y0, x1, y1)
            draw(tile, [(n, o.intersection(win)) for n, o in old if o.intersects(win)], 1.0, (x0, y0), OLD_COLOR, 3)
            draw(tile, [(n, o.intersection(win)) for n, o in new if o.intersects(win)], 1.0, (x0, y0), NEW_COLOR, 3)
            ImageDraw.Draw(tile).text((8, 8), f'{kind}: {area * PIXEL_UM ** 2:,.0f} um2 at ({cx * PIXEL_UM:.0f}, '
                                              f'{cy * PIXEL_UM:.0f}) um', fill=(255, 255, 255))
            tiles.append(tile)
        if tiles:
            w, h = max(t.width for t in tiles), max(t.height for t in tiles)
            grid = Image.new('RGB', (3 * w, 2 * h))
            for i, t in enumerate(tiles):
                grid.paste(t, ((i % 3) * w, (i // 3) * h))
            grid.save(OUTDIR / f'{tag}_tiles.png')


def one(tag: str) -> list[dict]:
    new, old = objects(NEW, tag), objects(OLD, tag)
    un, uo = unary_union([g for _, g in new]), unary_union([g for _, g in old])
    dice = 2 * un.intersection(uo).area / (un.area + uo.area)
    render(tag, new, old)
    return [{'image': tag, 'version': 'reviewed 2026-10-16 (run -20261015)', 'dice': dice, **metrics(tag, OLD, old)},
            {'image': tag, 'version': 'current 2026-10-27 (run -20261027)', 'dice': dice, **metrics(tag, NEW, new)}]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(a.workers) as ex:
        rows = [r for rs in ex.map(one, a.images) for r in rs]
    df = pd.DataFrame(rows)
    df.to_csv(OUTDIR / 'metrics.csv', index=False)
    pd.set_option('display.width', 250)
    print(df.round(3).to_string(index=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
