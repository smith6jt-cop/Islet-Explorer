#!/usr/bin/env python3
"""Review tiles for the v16 typing (owner standing rule: image-derived calls are checked against the image).

For an image and each class below: cells whose deciding marker stands just above its positivity level (0-0.5 widths)
and cells well above it (> 2 widths), three of each at random (seed 0), as 60 um tiles at full resolution: DAPI and the
deciding channels, each windowed from the measured background (black = the image's background for that marker:
the background peak, or the empty-slide null for E-cadherin and panCK, in image units) to the cohort scale's top
(p99.9), gamma 1; the cell is the circle at the center. For a top-marker label the "deciding" marker is the top marker,
and the two rows are cells above its level and cells named by it while below the level.

Reads the run's typing outputs (`typing/calls`, `typing/immune`, `typing/positivity_levels.csv`, `typing/shifts.json`,
`typing/scale.parquet`) and `cells/<image>.parquet` for centroids; no root is needed.
Writes output/typing_review_tiles_20261028/<image>/<class>.png.

    python scripts/research/typing_review_tiles_20261028.py --images ND_6232 T1D_6550
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from PIL import Image, ImageDraw, ImageFont

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform'))
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402
import redraw_tissue as RT  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'typing_review_tiles_20261028'
TILE_UM = 60.0
ZOOM = 3
PER_ROW = 3
# class -> (how to select, deciding marker, channels shown)
CLASSES = {
    'T Cell CD8a+': (('rule', 'Imm-T-cell', 'T-CD8'), 'CD3e', ['DAPI', 'CD3e', 'CD8a', 'CD4', 'CD45']),
    'T Cell CD4+': (('rule', 'Imm-T-cell', 'T-CD4-helper'), 'CD3e', ['DAPI', 'CD3e', 'CD4', 'CD8a', 'CD45']),
    'B Cell': (('rule', 'Imm-B-cell', None), 'CD20', ['DAPI', 'CD20', 'CD45', 'CD3e', 'ECAD']),
    'Macrophage': (('rule', 'Imm-macrophage', None), 'CD68', ['DAPI', 'CD68', 'CD163', 'CD45', 'ECAD']),
    'CD45+ (rule)': (('rule', 'Imm-CD45-leukocyte', None), 'CD45', ['DAPI', 'CD45', 'CD3e', 'CD68', 'ECAD']),
    'Exocrine (rule)': (('exocrine',), 'ECAD', ['DAPI', 'ECAD', 'panCK', 'CD45']),
    'SMA+ (top marker)': (('top', 'SMA+'), 'SMA', ['DAPI', 'SMA', 'ECAD', 'CD31']),
    'CD31+ (top marker)': (('top', 'CD31+'), 'CD31', ['DAPI', 'CD31', 'CD34', 'SMA']),
    'B3TUBB+ (top marker)': (('top', 'B3TUBB+'), 'B3TUBB', ['DAPI', 'B3TUBB', 'ECAD', 'CD45']),
    'INS+ (top marker)': (('top', 'INS+'), 'INS', ['DAPI', 'INS', 'CHGA', 'ECAD']),
    'CD45+ (top marker)': (('top', 'CD45+ (top marker)'), 'CD45', ['DAPI', 'CD45', 'ECAD', 'CD3e']),
}


def load_typing(tag: str) -> pd.DataFrame:
    typ = CT.typing_dir(RUN)
    calls = pd.read_parquet(typ / 'calls' / f'{tag}.parquet',
                            columns=['object_id', 'cell_class', 'immune_definition', 'resolved_by'])
    imm = pd.read_parquet(typ / 'immune' / f'{tag}.parquet', columns=['object_id', 'immune_class', 'immune_subtype'])
    xy = pd.read_parquet(RUN / 'cells' / f'{tag}.parquet', columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
    xy = xy.rename(columns={'Object ID': 'object_id', 'Centroid X µm': 'x', 'Centroid Y µm': 'y'})
    al = pd.read_parquet(typ / 'aligned' / f'{tag}.parquet')
    d = calls.merge(imm, on='object_id', how='left').merge(xy, on='object_id').merge(al, on='object_id')
    levels = CL.load_levels(RUN, tag)
    markers = [m for m in levels.index if m in d.columns]
    ex = CL.excess(d, markers, levels)
    ex.columns = [f'ex_{m}' for m in markers]
    return pd.concat([d.reset_index(drop=True), ex.reset_index(drop=True)], axis=1)


def select(d: pd.DataFrame, spec, marker: str, rng) -> tuple[pd.DataFrame, pd.DataFrame, str, str]:
    kind = spec[0]
    if kind == 'rule':
        use = d['immune_class'].fillna('').where(d['immune_class'].fillna('') != '', d['immune_definition'])
        m = (d['cell_class'] == 'Immune') & (use == spec[1]) & ~d['resolved_by'].str.endswith('top_marker')
        if spec[2]:
            m &= d['immune_subtype'] == spec[2]
    elif kind == 'exocrine':
        m = d['resolved_by'] == 'definition_exocrine'
    else:
        m = (d['cell_class'] == spec[1]) & d['resolved_by'].str.endswith('top_marker')
    g = d[m]
    e = g[f'ex_{marker}']
    if kind == 'top':
        a, b = g[e > 0], g[e <= 0]
        la, lb = f'{marker} above its level', f'{marker} below its level (still the top marker)'
    else:
        a, b = g[(e > 0) & (e <= 0.5)], g[e > 2]
        la, lb = f'{marker} 0-0.5 widths above its level', f'{marker} > 2 widths above its level'
    pick = lambda x: x.iloc[rng.choice(len(x), min(PER_ROW, len(x)), replace=False)] if len(x) else x
    return pick(a), pick(b), la, lb


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--images', nargs='+', required=True)
    a = ap.parse_args()
    typ = CT.typing_dir(RUN)
    scale = CL.load_scale(RUN)
    shifts = json.loads((typ / 'shifts.json').read_text())['shifts']
    font = ImageFont.load_default(size=16)
    half = int(round(TILE_UM / 2 / RT.PIXEL_UM))
    for tag in a.images:
        d = load_typing(tag)
        levels = CL.load_levels(RUN, tag)
        rng = np.random.default_rng(0)
        needed = sorted({c for _, _, chans in CLASSES.values() for c in chans})
        with tifffile.TiffFile(RT.IMAGES / f'{tag}.ome.tiff') as t:
            names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
            lev = t.series[0].levels[0]
            pages = {c: lev.pages[names.index(c)].asarray() for c in needed if c in names}
        H, W = next(iter(pages.values())).shape
        (OUT / tag).mkdir(parents=True, exist_ok=True)
        for cls, (spec, marker, chans) in CLASSES.items():
            near, far, l_near, l_far = select(d, spec, marker, rng)
            rows = [(l_near, near), (l_far, far)]
            n_rows = sum(len(x) for _, x in rows)
            if n_rows == 0:
                continue
            size = 2 * half * ZOOM
            canvas = Image.new('RGB', (len(chans) * (size + 6) + 6, n_rows * (size + 40) + 30), (20, 20, 20))
            draw = ImageDraw.Draw(canvas)
            for j, c in enumerate(chans):
                draw.text((6 + j * (size + 6), 6), c, fill=(255, 255, 255), font=font)
            r = 0
            for label, cells in rows:
                for _, cell in cells.iterrows():
                    cx, cy = int(round(cell['x'] / RT.PIXEL_UM)), int(round(cell['y'] / RT.PIXEL_UM))
                    y0, x0 = max(cy - half, 0), max(cx - half, 0)
                    y1, x1 = min(cy + half, H), min(cx + half, W)
                    for j, c in enumerate(chans):
                        if c not in pages:
                            continue
                        crop = pages[c][y0:y1, x0:x1].astype(np.float32)
                        f = shifts.get(c, {}).get('scale', {}).get(tag)
                        f = f if isinstance(f, float) else 1.0
                        if c == 'DAPI':
                            lo, hi = float(np.percentile(pages[c][::20, ::20], 1)), float(np.percentile(pages[c][::20, ::20], 99.7))
                        else:
                            lo = float(np.expm1(levels.loc[c, 'center'])) / f if c in levels.index and np.isfinite(levels.loc[c, 'center']) else 0.0
                            hi = float(scale.loc[c, 'top']) / f if c in scale.index else float(crop.max())
                        img = np.clip((crop - lo) / max(hi - lo, 1e-6), 0, 1)
                        tile = Image.fromarray((img * 255).astype(np.uint8)).resize((size, size), Image.NEAREST)
                        tile = tile.convert('RGB')
                        td = ImageDraw.Draw(tile)
                        rad = int(round(5.0 / RT.PIXEL_UM * ZOOM))
                        ctr = ((cx - x0) * ZOOM, (cy - y0) * ZOOM)
                        td.ellipse([ctr[0] - rad, ctr[1] - rad, ctr[0] + rad, ctr[1] + rad], outline=(255, 60, 60), width=2)
                        if c in levels.index and np.isfinite(cell.get(f'ex_{c}', np.nan)):
                            td.text((4, size - 22), f'{cell[f"ex_{c}"]:+.1f} w', fill=(255, 255, 0), font=font)
                        canvas.paste(tile, (6 + j * (size + 6), 30 + r * (size + 40)))
                    draw.text((6, 30 + r * (size + 40) + size + 4), f'{label}   x={cell["x"]:.0f} y={cell["y"]:.0f} um',
                              fill=(200, 200, 200), font=font)
                    r += 1
            canvas.save(OUT / tag / f'{cls.replace(" ", "_").replace("+", "pos").replace("(", "").replace(")", "")}.png')
        print(tag, 'done', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
