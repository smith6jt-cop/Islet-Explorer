#!/usr/bin/env python3
"""Review tiles for the final typing run (plan of 2026-10-28, Step 6): every label and status of the new hierarchy, on
every image, one sheet per label (the established endocrine labels are the owner-verified typing and are not re-shown).

A sheet has one row per image and two cells per row, a random one (seed 0) and the one whose naming marker stands
closest above its level (rows of four images). Each cell is a 48 um composite at full resolution, magnified 2x: DAPI gray, the naming marker
green, and magenta the cell's strongest positive marker of another compartment (or, if none, the marker that confirms
the label's lineage); each channel windowed from the tile's median to the image's p99.9, gamma 1. The cell is the white
circle; the caption gives its label, compartment status and the folds over level of the two colored markers.

Reads typing/identity/<image>.parquet, typing/context, the registry; writes
output/typing_final_checks_20261028/tiles_final/<label>.png.

    python scripts/research/typing_final_tiles_20261028.py --workers 8
"""
from __future__ import annotations

import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from PIL import Image, ImageDraw, ImageFont

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import registry as RG  # noqa: E402

RUN = REPO / 'data/development/uniform-20donor-20261027'
IMAGES = REPO / 'data/development/uniform-20donor-20261022-v10/images'
OUT = REPO / 'output/typing_final_checks_20261028/tiles_final'
PIXEL_UM = 0.3774
TILE_UM = 48.0
ZOOM = 2
PER_ROW = 4                      # images per sheet row, two cells each
STATUS_SHEETS = ('lineage_ambiguous', 'doublet_suspected')


def sheet_key(row) -> str:
    if row['named_by'] == 'established':
        return ''
    lab = str(row['label'])
    if 'Macrophage_' in lab:
        return 'Macrophage (all profiles)'
    return lab


def pick(tag: str) -> pd.DataFrame:
    """The cells shown for this image: per sheet a random cell and the cell closest above its naming marker's level."""
    d = pd.read_parquet(RUN / 'typing' / 'identity' / f'{tag}.parquet')
    xy = pd.read_parquet(RUN / 'typing' / 'context' / f'{tag}.parquet', columns=['object_id', 'x_um', 'y_um'])
    d = d.merge(xy, on='object_id', how='left')
    d['sheet'] = d.apply(sheet_key, axis=1)
    rng = np.random.default_rng(0)
    rows = []
    sheets = [(k, d['sheet'] == k) for k in d['sheet'].unique() if k] + \
             [(f'status {s}', d['compartment_status'] == s) for s in STATUS_SHEETS]
    for name, mask in sheets:
        g = d[mask]
        if not len(g):
            continue
        r = g.iloc[rng.integers(len(g))]
        rows.append({**r.to_dict(), 'sheet': name, 'which': 'random'})
        m = str(r['named_by']).split(';')[0]
        col = f'fold__{m}'
        if col in g.columns:
            near = g[g[col] > 0].nsmallest(1, col)
            if len(near):
                rows.append({**near.iloc[0].to_dict(), 'sheet': name, 'which': 'nearest its level'})
    return pd.DataFrame(rows)


def colors(reg: RG.Registry, r: pd.Series, markers: list[str]) -> tuple[str, str]:
    """(green, magenta) markers for one cell."""
    nb = [m for m in str(r['named_by']).split(';') if m in markers]
    folds = {m: r.get(f'fold__{m}', np.nan) for m in markers}
    green = nb[0] if nb else max((m for m in markers if r.get(f'pos__{m}', 0)), key=lambda m: folds[m], default='')
    comp = reg.table.loc[green, 'compartment'] if green in reg.table.index else ''
    other = [m for m in markers if r.get(f'pos__{m}', 0) and reg.role(m) == 'defining'
             and reg.table.loc[m, 'compartment'] != comp]
    if other:
        return green, max(other, key=lambda m: folds[m])
    conf = [m for m in markers if str(r.get('lineage', '')) in RG.split(reg.table.loc[m, 'confirms'])
            and m != green] if 'confirms' in reg.table.columns else []
    return green, (conf[0] if conf else '')


def render(tag: str) -> list[dict]:
    reg = RG.load()
    sel = pick(tag)
    if not len(sel):
        return []
    markers = [c[len('pos__'):] for c in sel.columns if c.startswith('pos__')]
    sel[['green', 'magenta']] = sel.apply(lambda r: pd.Series(colors(reg, r, markers)), axis=1)
    need = sorted({'DAPI'} | {reg.table.loc[m, 'channel'] for m in set(sel['green']) | set(sel['magenta']) if m})
    half = int(round(TILE_UM / PIXEL_UM / 2))
    font = ImageFont.load_default(size=13)
    out = []
    with tifffile.TiffFile(IMAGES / f'{tag}.ome.tiff') as t:
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        lev = t.series[0].levels[0]
        pages = {c: lev.pages[names.index(c)].asarray() for c in need if c in names}
    hi = {c: float(np.percentile(a[::16, ::16], 99.9)) for c, a in pages.items()}
    H, W = pages['DAPI'].shape
    for _, r in sel.iterrows():
        cx, cy = int(round(r['x_um'] / PIXEL_UM)), int(round(r['y_um'] / PIXEL_UM))
        y0, x0, y1, x1 = max(cy - half, 0), max(cx - half, 0), min(cy + half, H), min(cx + half, W)
        rgb = np.zeros((y1 - y0, x1 - x0, 3), np.float32)
        for m, col in ((None, (0.45, 0.45, 0.45)), (r['green'], (0, 1, 0)), (r['magenta'], (1, 0, 1))):
            c = 'DAPI' if m is None else (reg.table.loc[m, 'channel'] if m else '')
            if c not in pages:
                continue
            crop = pages[c][y0:y1, x0:x1].astype(np.float32)
            v = np.clip((crop - np.median(crop)) / max(hi[c] - np.median(crop), 1.0), 0, 1)
            for k in range(3):
                rgb[..., k] += v * col[k]
        img = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8)).resize(((x1 - x0) * ZOOM, (y1 - y0) * ZOOM))
        dr = ImageDraw.Draw(img)
        rad = 5.0 / PIXEL_UM * ZOOM
        c0 = ((cx - x0) * ZOOM, (cy - y0) * ZOOM)
        dr.ellipse([c0[0] - rad, c0[1] - rad, c0[0] + rad, c0[1] + rad], outline=(255, 255, 255), width=2)
        f = lambda m: f"{m} {r.get(f'fold__{m}', np.nan):+.1f}" if m else '-'
        cap = [f"{tag}  {r['which']}", f"{r['label']}  [{r['compartment_status']}]",
               f"green {f(r['green'])}  magenta {f(r['magenta'])}"]
        for j, line in enumerate(cap):
            dr.text((4, 4 + 15 * j), line, fill=(255, 255, 0), font=font)
        out.append({'sheet': r['sheet'], 'image': tag, 'which': r['which'], 'tile': np.asarray(img)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    tags = [r['image_name'] for r in RG.csv.DictReader(RG.COHORT.open())]
    with ProcessPoolExecutor(a.workers) as ex:
        tiles = [t for ts in ex.map(render, tags) for t in ts]
    OUT.mkdir(parents=True, exist_ok=True)
    size = int(round(TILE_UM / PIXEL_UM / 2)) * 2 * ZOOM
    for sheet in sorted({t['sheet'] for t in tiles}):
        ts = [t for t in tiles if t['sheet'] == sheet]
        rows = [tag for tag in tags if any(t['image'] == tag for t in ts)]
        n_rows = -(-len(rows) // PER_ROW)
        canvas = Image.new('RGB', (2 * PER_ROW * (size + 4) + 4, n_rows * (size + 4) + 4), (25, 25, 25))
        for i, tag in enumerate(rows):
            for t in ts:
                if t['image'] != tag:
                    continue
                col = (i % PER_ROW) * 2 + (0 if t['which'] == 'random' else 1)
                canvas.paste(Image.fromarray(t['tile']), (4 + col * (size + 4), 4 + (i // PER_ROW) * (size + 4)))
        safe = re.sub(r'[^A-Za-z0-9_.+-]+', '_', sheet).replace('+', 'pos')
        canvas.save(OUT / f'{safe}.png')
    print(f'{len({t["sheet"] for t in tiles})} sheets in {OUT}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
