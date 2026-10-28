#!/usr/bin/env python3
"""Evidence for the final typing run's positivity step (plan of 2026-10-28, Step 2): histograms and image tiles.

- histograms/<marker>.png: for every marker with a computed level, the 20 images' histograms of log1p(REDSEA value) of
  analysis-eligible cells (exact zeros in the first bar), the chosen level (red) and the background peak (gray dashed),
  from typing/levels.csv; the marker's decision (use / drop / stop) in the title.
- tiles/CD20_<image>.png: CD20-positive cells at the chosen level, three CD45-negative and three CD45-positive, 60 um
  tiles at full resolution: DAPI dim gray, CD20 green, CD45 magenta, E-cadherin blue.
- tiles/SMA_<image>.png: SMA-positive endocrine-object hormone-positive cells and leukocytes with no other SMA-positive
  cell within 20 um, six per image: DAPI dim gray, SMA green, hormones (INS+GCG+SST) magenta, CD45 blue.
Each channel is windowed from the tile's tissue median (black) to the image's p99.9 of that channel, gamma 1; the cell
is the circle at the center.

    python scripts/research/typing_final_evidence_20261028.py --tiles-images ND_6548 ND_6356 T1D_6533 Aab_6538
"""
from __future__ import annotations

import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tifffile  # noqa: E402
from PIL import Image, ImageDraw  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import registry as RG  # noqa: E402
import levels as LV  # noqa: E402

RUN = REPO / 'data/development/uniform-20donor-20261027'
IMAGES = REPO / 'data/development/uniform-20donor-20261022-v10/images'
OUT = REPO / 'output/typing_final_checks_20261028'
PIXEL_UM = 0.3774
TILE_UM = 60.0
CHANNEL = {'PGP9.5': 'PGP9_5', 'KRT8-18': 'KRT8_18'}


def histograms(reg: RG.Registry) -> None:
    lv = pd.read_csv(RUN / 'typing' / 'levels.csv')
    sel = pd.read_csv(RUN / 'typing' / 'levels_selection.csv').set_index('marker')
    tags = sorted(lv['image'].unique())
    data = {t: LV.load(RUN, t, reg) for t in tags}
    (OUT / 'histograms').mkdir(parents=True, exist_ok=True)
    for m in sorted(lv['marker'].unique()):
        fig, axes = plt.subplots(4, 5, figsize=(20, 13), sharex=True)
        rows = lv[lv['marker'] == m].set_index('image')
        for ax, t in zip(axes.ravel(), tags):
            r = rows.loc[t]
            x = data[t][reg.column(m, r['measurement'])].to_numpy(float)
            v = np.log1p(np.clip(np.nan_to_num(x, nan=0.0), 0, None))
            ax.hist(v, bins=LV.EDGES, color='0.35')
            if np.isfinite(r['center']):
                ax.axvline(r['center'], color='0.6', ls='--', lw=1)
            if np.isfinite(r['level_log']):
                ax.axvline(r['level_log'], color='red', lw=1.5)
            ax.set_title(f"{t}: {100 * r['positive_fraction']:.1f} % positive", fontsize=10)
            ax.set_yscale('log')
        d = sel.loc[m]
        fig.suptitle(f"{m} ({d['measurement']}, {d['method']}): decision {d['decision']}; log1p(REDSEA value); "
                     f"red = level, dashed = background peak", fontsize=13)
        fig.tight_layout()
        fig.savefig(OUT / 'histograms' / f'{m}.png', dpi=110)
        plt.close(fig)


def window(a: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.clip((a.astype(np.float32) - lo) / max(hi - lo, 1.0), 0, 1)


def tiles(tag: str) -> None:
    reg = RG.load()
    lv = pd.read_csv(RUN / 'typing' / 'levels.csv')
    lv = lv[lv['image'] == tag].set_index('marker')
    d = LV.load(RUN, tag, reg)

    def pos(m):
        r = lv.loc[m]
        x = d[reg.column(m, r['measurement'])].to_numpy(float)
        return np.log1p(np.clip(np.nan_to_num(x, nan=0.0), 0, None)) > r['level_log']

    p = {m: pos(m) for m in ('CD20', 'CD45', 'SMA', 'INS', 'GCG', 'SST', 'CD3e', 'CD68', 'CD163')}
    xy = d[['x_um', 'y_um']].to_numpy(float)
    rng = np.random.default_rng(0)
    pick = lambda mask, n: rng.choice(np.flatnonzero(mask), min(n, int(mask.sum())), replace=False) if mask.any() else []
    leuk = p['CD45'] & (p['CD3e'] | p['CD68'] | p['CD163'])
    hormone = (p['INS'] | p['GCG'] | p['SST']) & d['in_endocrine'].to_numpy()
    src = np.flatnonzero(p['SMA'])
    nn = cKDTree(xy[src]).query(xy, k=2)[0][:, 1] if len(src) > 1 else np.full(len(d), np.inf)
    sets = {
        'CD20': ([pick(p['CD20'] & ~p['CD45'], 3), pick(p['CD20'] & p['CD45'], 3)],
                 {'DAPI': (0.3, 0.3, 0.3), 'CD20': (0, 1, 0), 'CD45': (1, 0, 0.85), 'ECAD': (0.25, 0.55, 1)},
                 ['CD20+ CD45-'] * 3 + ['CD20+ CD45+'] * 3),
        'SMA': ([pick(p['SMA'] & hormone & (nn > 20), 3), pick(p['SMA'] & leuk & (nn > 20), 3)],
                {'DAPI': (0.3, 0.3, 0.3), 'SMA': (0, 1, 0), 'INS': (1, 0, 0.85), 'GCG': (1, 0, 0.85), 'SST': (1, 0, 0.85),
                 'CD45': (0.25, 0.55, 1)},
                ['SMA+ endocrine cell, no SMA+ cell within 20 um'] * 3 + ['SMA+ leukocyte, no SMA+ cell within 20 um'] * 3),
    }
    with tifffile.TiffFile(IMAGES / f'{tag}.ome.tiff') as t:
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        lev = t.series[0].levels[0]
        H, W = lev.shape[-2:]
        half = int(TILE_UM / PIXEL_UM / 2)
        cache = {}
        for name, (groups, colors, captions) in sets.items():
            idx = [i for g in groups for i in g]
            if not idx:
                continue
            panels = []
            for i, cap in zip(idx, captions[:len(idx)] if len(idx) == 6 else captions):
                cx, cy = xy[i] / PIXEL_UM
                x0, y0, x1, y1 = int(max(0, cx - half)), int(max(0, cy - half)), int(min(W, cx + half)), int(min(H, cy + half))
                rgb = np.zeros((y1 - y0, x1 - x0, 3), np.float32)
                for ch, col in colors.items():
                    if ch not in cache:
                        a = lev.pages[names.index(CHANNEL.get(ch, ch))].asarray()
                        cache[ch] = (a, float(np.percentile(a[::8, ::8], 99.9)))
                    a, hi = cache[ch]
                    crop = a[y0:y1, x0:x1]
                    w = window(crop, float(np.median(crop)), hi)
                    for k in range(3):
                        rgb[..., k] += w * col[k]
                img = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8)).resize(((x1 - x0) * 2, (y1 - y0) * 2))
                dr = ImageDraw.Draw(img)
                r = 2 * 6 / PIXEL_UM
                c = (2 * (cx - x0), 2 * (cy - y0))
                dr.ellipse([c[0] - r, c[1] - r, c[0] + r, c[1] + r], outline=(255, 255, 255), width=2)
                dr.text((6, 6), cap, fill=(255, 255, 255))
                panels.append(img)
            w_, h_ = max(i.width for i in panels), max(i.height for i in panels)
            grid = Image.new('RGB', (3 * w_, 2 * h_))
            for j, im in enumerate(panels):
                grid.paste(im, ((j % 3) * w_, (j // 3) * h_))
            (OUT / 'tiles').mkdir(parents=True, exist_ok=True)
            grid.save(OUT / 'tiles' / f'{name}_{tag}.png')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tiles-images', nargs='*', default=[])
    ap.add_argument('--no-histograms', action='store_true')
    a = ap.parse_args()
    if not a.no_histograms:
        histograms(RG.load())
    if a.tiles_images:
        with ProcessPoolExecutor(len(a.tiles_images)) as ex:
            list(ex.map(tiles, a.tiles_images))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
