#!/usr/bin/env python3
"""Background from the empty (unannotated, unsegmented) slide area, per image and channel (owner question,
2026-10-28: "have you tried using the background from the empty areas of the image that were not annotated or
segmented?"). Measurements only.

Per image:
  - empty area: full-resolution pixels at least 20 um outside the image's Tissue outline (the run's
    `inputs/tissue/<image>.geojson`), inside the image;
  - the null for a cell: the mean of a WINDOW x WINDOW-pixel square of empty area (about one cell's measurement area),
    sampled N_WINDOWS times at random positions wholly inside the empty area; per channel its p50 / p99 / p99.9;
  - cells: each marker's REDSEA-corrected value in its ruled compartment (`typing/redsea/<image>.parquet`,
    before UniFORM) and QuPath's uncorrected value in the same compartment (`cells/<image>.parquet`);
  - positivity = above the empty null's p99.9; the label-free consistency checks of
    `typing_accuracy_checks_20261028.py` on those calls.
Written to output/empty_area_background_20261028/.

    python scripts/research/empty_area_background_20261028.py
"""
from __future__ import annotations

import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
from shapely.geometry import shape

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform'))
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import redraw_tissue as RT  # noqa: E402
import codex_typing as CT  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'empty_area_background_20261028'
MARKERS = ['CD3e', 'CD8a', 'CD4', 'CD20', 'CD45', 'CD68', 'CD163', 'ECAD', 'panCK', 'BCatenin', 'SMA', 'B3TUBB', 'CD31', 'VIM']
WINDOW = 19            # px at 0.3774 um: 7.2 um square, about a cell's measurement area
N_WINDOWS = 20000
MARGIN_UM = 20.0
SEED = 0


def empty_null(tag: str) -> pd.DataFrame:
    raw = json.loads((RUN / 'inputs' / 'tissue' / f'{tag}.geojson').read_text())
    geom = shape((raw['features'] if isinstance(raw, dict) else raw)[0]['geometry'])
    with tifffile.TiffFile(RT.IMAGES / f'{tag}.ome.tiff') as t:
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        lev = t.series[0].levels[0]
        H, W = lev.shape[-2:]
        mask = np.ones((H, W), np.uint8)
        for p in getattr(geom, 'geoms', [geom]):
            cv2.fillPoly(mask, [np.round(np.asarray(p.exterior.coords)).astype(np.int32)], 0)
            for h in p.interiors:
                cv2.fillPoly(mask, [np.round(np.asarray(h.coords)).astype(np.int32)], 1)
        r = int(round(MARGIN_UM / RT.PIXEL_UM))
        mask = cv2.erode(mask, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1)))
        # window centers whose whole window is empty area
        ok = cv2.erode(mask, np.ones((WINDOW, WINDOW), np.uint8)).astype(bool)
        e = WINDOW // 2 + 1                                         # the whole window inside the image
        ok[:e], ok[-e:], ok[:, :e], ok[:, -e:] = False, False, False, False
        ys, xs = np.nonzero(ok)
        rows = []
        if len(ys) == 0:
            return pd.DataFrame([{'image': tag, 'marker': m, 'windows': 0} for m in MARKERS])
        rng = np.random.default_rng(SEED)
        pick = rng.choice(len(ys), min(N_WINDOWS, len(ys)), replace=False)
        ys, xs = ys[pick], xs[pick]
        h = WINDOW // 2
        for m in MARKERS:
            a = lev.pages[names.index(m)].asarray().astype(np.float64)
            ii = cv2.integral(a)
            y0, y1, x0, x1 = ys - h, ys + h + 1, xs - h, xs + h + 1
            s = ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]
            means = s / WINDOW ** 2
            rows.append({'image': tag, 'marker': m, 'windows': len(means), 'empty_area_px': int(mask.sum()),
                         'null_p50': np.percentile(means, 50), 'null_p99': np.percentile(means, 99),
                         'null_p999': np.percentile(means, 99.9), 'pixel_p50': float(np.median(a[mask.astype(bool)]))})
    return pd.DataFrame(rows)


def cell_values(tag: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    typ = CT.typing_dir(RUN)
    rs = pd.read_parquet(typ / 'redsea' / f'{tag}.parquet')
    rs = rs[rs['qc_analysis_eligible']]
    red = pd.DataFrame({m: rs[f'{CT.COMPARTMENT[m]}__{m}'].to_numpy(float) for m in MARKERS})
    for m in MARKERS:      # owner, 2026-10-27: no cell without cytoplasm (as align_scale.read_measured)
        if CT.COMPARTMENT[m] == 'Cytoplasm':
            red[m] = red[m].where(red[m].notna(), rs[f'Cell__{m}'].to_numpy(float))
    qp = pd.read_parquet(RUN / 'cells' / f'{tag}.parquet')
    qp = qp.set_index('Object ID').reindex(rs['object_id'])
    raw = pd.DataFrame({m: qp[f'{CT.COMPARTMENT[m]}: {m}: Mean'].to_numpy(float)
                        if f'{CT.COMPARTMENT[m]}: {m}: Mean' in qp else np.nan for m in MARKERS})
    return red, raw


def checks(P: pd.DataFrame) -> dict:
    c3 = P['CD3e']
    return {**{f'pos%_{m}': 100 * P[m].mean() for m in MARKERS},
            'CD45+|CD3e+': P['CD45'][c3].mean(), 'CD3e+|CD8a+': c3[P['CD8a']].mean(),
            'CD20+|CD3e+': P['CD20'][c3].mean(), 'CD68+|CD3e+': P['CD68'][c3].mean(),
            'CD45+|CD20+': P['CD45'][P['CD20']].mean(), 'CD45+|CD68+': P['CD45'][P['CD68']].mean()}


def one(tag: str) -> tuple[pd.DataFrame, list[dict]]:
    null = empty_null(tag)
    red, raw = cell_values(tag)
    nl = null.set_index('marker')
    out = []
    for name, vals in (('REDSEA', red), ('QuPath uncorrected', raw)):
        for q in ('null_p99', 'null_p999'):
            if 'null_p999' not in nl or nl['windows'].min() == 0:
                continue
            P = vals.gt(nl.loc[MARKERS, q].to_numpy(), axis=1)
            out.append({'image': tag, 'values': name, 'level': q, **checks(P)})
    print(tag, 'windows', int(null['windows'].min()), flush=True)
    return null, out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    tags = [r['image_name'] for r in CT.cohort_rows()]
    nulls, recs = [], []
    with ProcessPoolExecutor(5) as ex:
        for null, rec in ex.map(one, tags):
            nulls.append(null)
            recs += rec
    pd.concat(nulls).to_csv(OUT / 'empty_area_null_by_image.csv', index=False)
    R = pd.DataFrame(recs)
    R.to_csv(OUT / 'checks_by_image.csv', index=False)
    pd.set_option('display.width', 250)
    print(R.groupby(['values', 'level']).median(numeric_only=True).round(3).T.to_string())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
