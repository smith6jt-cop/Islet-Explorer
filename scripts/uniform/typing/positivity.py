#!/usr/bin/env python3
"""positivity.py -- each marker's positivity level per image: its own background on the log scale.

Owner, 2026-10-28 ("almost every other pipeline uses normalization and scaling"; approved plan): positivity is decided
on log1p(value) per image, as the published CODEX / MIBI pipelines do (a log or arcsinh transform, then per-sample
standardization), not at a fraction of a linear cohort range: v15's level, scaled >= 0.10, sat at 3-170x the image
background depending on the marker, and 2.6-19x across images for one marker
(output/positivity_method_comparison_20261028/).

Markers most cells lack (the background peak is the negative population):
  center    the mode of log1p(x) over the image's analysis-eligible cells with x > 0 (x is the aligned own-compartment
            REDSEA value; UniFORM's factor is an additive shift on the log scale that the per-image center absorbs)
  width     the peak's left half-width at half-maximum / 1.1774 (a Gaussian's sigma), pooled over the 20 images as
            the median per marker (codex_typing.POSITIVITY_WIDTH_POOLED): where REDSEA's clamp at 0 removed 25-46 %
            of the cells the remaining peak is a plateau and its own width is unreadable (levels 18-196x
            background), and the widths agree where readable (CD68 0.47-0.71 log units)
  level     center + POSITIVITY_K x width (K = 3: three background widths above the peak, the usual convention)
Markers most cells carry (codex_typing.MAJORITY_MARKERS; no background peak in tissue):
  center    the empty-slide null: the EMPTY_QUANTILE percentile of the mean of a cell-sized window of image pixels
            at least EMPTY_MARGIN_UM outside the Tissue outline, times the image's UniFORM factor (owner: "the
            background outside the tissue [helps] with the high expressers like ECAD")
  level     from the marker's own expressing population, floored at the null (codex_typing.MAJORITY_LEVEL):
            'edge' = the population's peak minus 3 right-side widths (the side with no non-epithelial cell;
            E-cadherin), 'gate' = the two-component mixture gate on the nonzero log values (bimodal; panCK)
  width     the expressing population's right-side width (1.0 where unreadable)
excess (classify.excess): (log1p(x) - level) / width, > 0 positive, in widths above the level; it also orders the
top-marker step's markers.

Written to `$OUT/typing/positivity_levels.csv` (image, marker, center, width, image_width, level_log, level_raw,
source, n, zero_fraction). Run after align_scale.py; classify.py, immune_subtype.py and assemble_uniform.py read it.

    python scripts/uniform/typing/positivity.py --out $OUT
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
from shapely.geometry import shape
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402
import redraw_tissue as RT  # noqa: E402


def peak(v: np.ndarray, side: str = 'left') -> tuple[float, float]:
    """The peak of log values: its mode, and the sigma from its half-width at half-maximum on one side (5-bin
    smoothed 256-bin histogram); sigma is NaN when no bin on that side falls under half the peak's height."""
    v = v[np.isfinite(v)]
    if len(v) < 1000:
        return np.nan, np.nan
    h, e = np.histogram(v, bins=256)
    hs = np.convolve(h.astype(float), np.ones(5) / 5, mode='same')
    i = int(np.argmax(hs[1:])) + 1
    c = (e[:-1] + e[1:]) / 2
    if side == 'left':
        below = np.flatnonzero(hs[:i] <= hs[i] / 2)
        return float(c[i]), (float((c[i] - c[int(below[-1])]) / 1.1774) if len(below) else np.nan)
    above = np.flatnonzero(hs[i:] <= hs[i] / 2)
    return float(c[i]), (float((c[i + int(above[0])] - c[i]) / 1.1774) if len(above) else np.nan)


def mixture_gate(v: np.ndarray, n: int = 200000, seed: int = CT.EMPTY_SEED) -> float:
    """scimap.pp.rescale's gmm gate: GaussianMixture(2) on the log values, the mean of the two component means."""
    v = v[np.isfinite(v)]
    if len(v) < 1000:
        return np.nan
    if len(v) > n:
        v = np.random.default_rng(seed).choice(v, n, replace=False)
    return float(np.mean(GaussianMixture(n_components=2, random_state=seed).fit(v.reshape(-1, 1)).means_))


def image_peaks(out: Path, tag: str, markers: list[str]) -> tuple[str, dict]:
    al = pd.read_parquet(CT.typing_dir(out) / 'aligned' / f'{tag}.parquet')
    al = al[al['qc_analysis_eligible']]
    rows = {}
    for m in markers:
        if m not in al.columns:                 # not acquired on this image (codex_typing.NOT_ON_EVERY_IMAGE)
            rows[m] = {'center': np.nan, 'width': np.nan, 'n': 0, 'zeros': np.nan}
            continue
        x = al[m].to_numpy(float)
        v = np.log1p(np.clip(x, 0, None))[x > 0]
        if m in CT.MAJORITY_MARKERS:
            mode, right = peak(v, 'right')
            rows[m] = {'mode': mode, 'width': right, 'gate': mixture_gate(v), 'n': int(np.isfinite(x).sum()),
                       'zeros': float(np.mean(x <= 0))}
        else:
            center, width = peak(v, 'left')
            rows[m] = {'center': center, 'width': width, 'n': int(np.isfinite(x).sum()), 'zeros': float(np.mean(x <= 0))}
    return tag, rows


def empty_null(out: Path, tag: str, markers: list[str]) -> tuple[str, dict]:
    """The EMPTY_QUANTILE percentile of cell-sized windows of empty slide (outside the Tissue outline), per channel,
    in the image's own units."""
    raw = json.loads((out / 'inputs' / 'tissue' / f'{tag}.geojson').read_text())
    geom = shape((raw['features'] if isinstance(raw, dict) else raw)[0]['geometry'])
    with tifffile.TiffFile(RT.IMAGES / f'{tag}.ome.tiff') as t:
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        lev = t.series[0].levels[0]
        H, W = lev.shape[-2:]
        mask = np.ones((H, W), np.uint8)
        for p in getattr(geom, 'geoms', [geom]):
            cv2.fillPoly(mask, [np.round(np.asarray(p.exterior.coords)).astype(np.int32)], 0)
            for hole in p.interiors:
                cv2.fillPoly(mask, [np.round(np.asarray(hole.coords)).astype(np.int32)], 1)
        r = int(round(CT.EMPTY_MARGIN_UM / CT.PIXEL_UM))
        mask = cv2.erode(mask, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1)))
        w = CT.EMPTY_WINDOW_PX
        ok = cv2.erode(mask, np.ones((w, w), np.uint8)).astype(bool)
        e = w // 2 + 1
        ok[:e], ok[-e:], ok[:, :e], ok[:, -e:] = False, False, False, False
        ys, xs = np.nonzero(ok)
        if len(ys) < 100:
            return tag, {m: np.nan for m in markers}
        rng = np.random.default_rng(CT.EMPTY_SEED)
        pick = rng.choice(len(ys), min(CT.EMPTY_WINDOWS, len(ys)), replace=False)
        ys, xs = ys[pick], xs[pick]
        h = w // 2
        null = {}
        for m in markers:
            if m not in names:
                null[m] = np.nan
                continue
            a = lev.pages[names.index(m)].asarray().astype(np.float64)
            ii = cv2.integral(a)
            s = ii[ys + h + 1, xs + h + 1] - ii[ys - h, xs + h + 1] - ii[ys + h + 1, xs - h] + ii[ys - h, xs - h]
            null[m] = float(np.percentile(s / w ** 2, CT.EMPTY_QUANTILE))
    return tag, null


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=5)
    a = ap.parse_args()
    out = a.out.resolve()
    typ = CT.typing_dir(out)
    scale = CL.load_scale(out)
    shifts = json.loads((typ / 'shifts.json').read_text())['shifts']
    tags = [r['image_name'] for r in CT.cohort_rows()]
    markers = list(scale.index)
    majority = [m for m in markers if m in CT.MAJORITY_MARKERS]
    with ProcessPoolExecutor(a.workers) as ex:
        peaks = dict(ex.map(image_peaks, repeat(out), tags, repeat(markers)))
        nulls = dict(ex.map(empty_null, repeat(out), tags, repeat(majority)))
    pooled = {m: float(np.nanmedian([peaks[t][m]['width'] for t in tags])) for m in markers if m not in majority}
    rows = []
    for t in tags:
        for m in markers:
            p = peaks[t][m]
            if m in majority:
                f = shifts.get(m, {}).get('scale', {}).get(t)
                f = f if isinstance(f, float) else 1.0
                null = np.log1p(nulls[t][m] * f) if np.isfinite(nulls[t].get(m, np.nan)) else np.nan
                how = CT.MAJORITY_LEVEL.get(m, 'edge')
                own = p.get('gate', np.nan) if how == 'gate' else p.get('mode', np.nan) - CT.POSITIVITY_K * p.get('width', np.nan)
                level = np.nanmax([own, null]) if np.isfinite(own) or np.isfinite(null) else np.nan
                width = p.get('width', np.nan) if np.isfinite(p.get('width', np.nan)) else 1.0
                source = ('mixture_gate' if how == 'gate' else 'epithelial_edge') if np.isfinite(own) else \
                    ('empty_area' if np.isfinite(null) else 'unavailable')
                rows.append({'image': t, 'marker': m, 'center': null, 'width': width, 'image_width': p.get('width', np.nan),
                             'level_log': level, 'level_raw': float(np.expm1(level)) if np.isfinite(level) else np.nan,
                             'source': source, 'n': p['n'], 'zero_fraction': p['zeros']})
            else:
                w = pooled[m] if CT.POSITIVITY_WIDTH_POOLED else p['width']
                level = p['center'] + CT.POSITIVITY_K * w if np.isfinite(p['center']) and np.isfinite(w) else np.nan
                rows.append({'image': t, 'marker': m, 'center': p['center'], 'width': w, 'image_width': p['width'],
                             'level_log': level, 'level_raw': float(np.expm1(level)) if np.isfinite(level) else np.nan,
                             'source': 'background_peak' if np.isfinite(level) else 'no_peak', 'n': p['n'],
                             'zero_fraction': p['zeros']})
    df = pd.DataFrame(rows)
    df.to_csv(typ / 'positivity_levels.csv', index=False)
    fold = (df['level_raw'] / np.expm1(df['center'])).where(df['source'] == 'background_peak')
    summary = df.assign(fold_over_background=fold).groupby('marker').agg(
        source=('source', 'first'), width=('width', 'first'), level_raw_median=('level_raw', 'median'),
        fold_over_background=('fold_over_background', 'median'), zero_fraction_max=('zero_fraction', 'max'))
    print(summary.round(3).to_string())
    print('sources:', df['source'].value_counts().to_dict(), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
