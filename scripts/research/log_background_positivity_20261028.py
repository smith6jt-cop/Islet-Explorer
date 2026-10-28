#!/usr/bin/env python3
"""Positivity on log-transformed, per-image standardized values (owner, 2026-10-28: "almost every other pipeline uses
normalization and scaling"). Measurements only.

Per image and marker, on the REDSEA own-compartment values (whole-cell fallback where the nucleus fills the cell):
  v = log1p(x)                                      the multiplicative scale the published pipelines use
  background center = the mode of v over x > 0      (the most common value, as align_scale.py finds the floor)
and a background spread, three estimates:
  lower half        SD of v below the mode, mirrored. Inflated by REDSEA's subtraction tail (values pushed toward 0
                    sit far below the mode on the log scale): levels 20-50x background. Kept for the record.
  peak width        the background peak's left half-width at half-maximum / 1.1774 (a Gaussian's sigma), reading
                    only the peak. Readable where few cells were clamped to 0; where REDSEA's clamp removed 25-46 %
                    of the cells the remaining peak is distorted and the level lands 18-196x background.
  cohort width      the median peak width over the 20 images, per marker, applied to each image's own center: a
                    per-image background with a pooled dispersion, the usual remedy for an unreadable per-sample
                    width.
A cell is positive at z = (v - center) / spread > K (K = 3, 4, 5). Also:
  scimap gate       the method the original 15-donor typing used on this dataset (scimap.pp.rescale, method gmm):
                    a two-component Gaussian mixture on log1p values per image and marker, gate = the mean of the
                    two component means, positive above the gate. Fitted on the nonzero values: with the clamped
                    zeros in, the mixture splits zero from everything else (70-90 % positive); without them it
                    still splits the background (20-58 % positive) except on bimodal markers such as CD31.
All are compared with v15's level (scaled >= 0.10 on the linear cohort scale) on the label-free checks of
typing_accuracy_checks_20261028.py. Markers most cells carry (E-cadherin, panCK) have no background peak in tissue
and are shown only for completeness: their floor is the empty-area null (empty_area_background_20261028.py).
A figure per image (figures/<image>_log_histograms.png): each marker's log1p histogram with its background peak,
the peak-width and cohort-width levels (z = 3), the scimap gate and v15's level.
Written to output/positivity_method_comparison_20261028/log_background_{checks,levels}.csv.

    python scripts/research/log_background_positivity_20261028.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'positivity_method_comparison_20261028'
MARKERS = ['CD3e', 'CD8a', 'CD4', 'CD20', 'CD45', 'CD68', 'CD163', 'LGALS3', 'SMA', 'B3TUBB', 'CD31', 'ECAD', 'panCK']
KS = (3, 4, 5)


def load(tag: str) -> pd.DataFrame:
    typ = CT.typing_dir(RUN)
    rs = pd.read_parquet(typ / 'redsea' / f'{tag}.parquet')
    rs = rs[rs['qc_analysis_eligible']]
    X = pd.DataFrame({m: rs[f'{CT.COMPARTMENT[m]}__{m}'].to_numpy(float) for m in MARKERS})
    for m in MARKERS:                      # owner, 2026-10-27: no cell without cytoplasm
        if CT.COMPARTMENT[m] == 'Cytoplasm':
            X[m] = X[m].where(X[m].notna(), rs[f'Cell__{m}'].to_numpy(float))
    X.insert(0, 'object_id', rs['object_id'].to_numpy())
    return X.reset_index(drop=True)


def background(v: np.ndarray) -> tuple[float, float]:
    """Mode and lower-half SD of log values (x > 0 only, as align_scale.py)."""
    v = v[np.isfinite(v)]
    h, e = np.histogram(v, bins=256)
    i = int(np.argmax(h[1:])) + 1
    mode = float((e[i] + e[i + 1]) / 2)
    low = v[v < mode]
    spread = float(np.sqrt(np.mean((mode - low) ** 2))) if len(low) > 100 else np.nan
    return mode, spread


def peak_width(v: np.ndarray) -> tuple[float, float]:
    """Mode and sigma from the background peak's left half-width at half-maximum (5-bin smoothed histogram)."""
    v = v[np.isfinite(v)]
    h, e = np.histogram(v, bins=256)
    hs = np.convolve(h.astype(float), np.ones(5) / 5, mode='same')
    i = int(np.argmax(hs[1:])) + 1
    c = (e[:-1] + e[1:]) / 2
    mode = float(c[i])
    below = np.flatnonzero(hs[:i] <= hs[i] / 2)
    if len(below) == 0:
        return mode, np.nan
    j = int(below[-1])
    return mode, float((c[i] - c[j]) / 1.1774)


def scimap_gate(v: np.ndarray, n: int = 200000, seed: int = 0) -> float:
    """scimap.pp.rescale's gmm gate: GaussianMixture(2) on log1p values, the mean of the two component means."""
    v = v[np.isfinite(v)]
    if len(v) > n:
        v = np.random.default_rng(seed).choice(v, n, replace=False)
    gmm = GaussianMixture(n_components=2, random_state=seed).fit(v.reshape(-1, 1))
    return float(np.mean(gmm.means_))


def figure(t: str, X: pd.DataFrame, lv: dict, path: Path) -> None:
    fig, axes = plt.subplots(3, 5, figsize=(24, 13))
    for ax, m in zip(axes.flat, MARKERS):
        x = X[m].to_numpy(float)
        v = np.log1p(np.clip(x, 0, None))
        ax.hist(v[x > 0], bins=200, color='#999999', density=True)
        r = lv[m]
        ax.axvline(r['peak_mode'], color='#4477AA', ls=':', lw=2, label='background peak')
        ax.axvline(np.log1p(r['peak_level_k3']), color='#228833', ls='--', lw=2, label='this image\'s peak width, z = 3')
        ax.axvline(np.log1p(r['cohort_level_k3']), color='#EE6677', lw=2, label='cohort width, z = 3')
        ax.axvline(np.log1p(r['scimap_gate']), color='#AA3377', ls='-.', lw=2, label='scimap gate (nonzero values)')
        ax.axvline(np.log1p(r['v15_level_linear']), color='black', lw=2, label='v15: 0.10 on the linear scale')
        ax.set_yscale('log')
        ax.set_title(f'{m}  ({100 * (x <= 0).mean():.0f} % of cells at 0)', fontsize=14)
        ax.set_xlabel('log1p(REDSEA value, own compartment)', fontsize=11)
    for ax in list(axes.flat)[len(MARKERS):]:
        ax.axis('off')
    axes.flat[0].legend(fontsize=10, loc='upper right')
    fig.suptitle(f'{t}: per-marker distributions (nonzero cells, log y) and the positivity levels compared', fontsize=17)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def checks(P: pd.DataFrame) -> dict:
    c3 = P['CD3e']
    return {**{f'pos%_{m}': 100 * P[m].mean() for m in MARKERS},
            'CD45+|CD3e+': P['CD45'][c3].mean(), 'CD3e+|CD8a+': c3[P['CD8a']].mean(),
            'CD20+|CD3e+': P['CD20'][c3].mean(), 'CD68+|CD3e+': P['CD68'][c3].mean(),
            'CD45+|CD20+': P['CD45'][P['CD20']].mean(), 'CD45+|CD68+': P['CD45'][P['CD68']].mean()}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'figures').mkdir(exist_ok=True)
    typ = CT.typing_dir(RUN)
    scale = CL.load_scale(RUN)
    tags = [r['image_name'] for r in CT.cohort_rows()]
    data = {t: load(t) for t in tags}
    # pass 1: each image's background peak and its width, per marker; the cohort width is their median
    peaks = {}
    for t in tags:
        for m in MARKERS:
            x = data[t][m].to_numpy(float)
            v = np.log1p(np.clip(x, 0, None))
            peaks[(t, m)] = peak_width(v[x > 0])
    cohort_spread = {m: float(np.nanmedian([peaks[(t, m)][1] for t in tags])) for m in MARKERS}
    # pass 2: levels, calls, checks, figures
    recs, levels = [], []
    for t in tags:
        X = data[t]
        Z = pd.DataFrame(index=X.index)        # lower-half spread
        Zp = pd.DataFrame(index=X.index)       # this image's peak width
        Zc = pd.DataFrame(index=X.index)       # cohort width
        G = pd.DataFrame(index=X.index)        # scimap gate
        for m in MARKERS:
            x = X[m].to_numpy(float)
            v = np.log1p(np.clip(x, 0, None))
            mode, spread = background(v[x > 0])
            mode_p, spread_p = peaks[(t, m)]
            gate = scimap_gate(v[x > 0])
            Z[m] = (v - mode) / spread
            Zp[m] = (v - mode_p) / spread_p
            Zc[m] = (v - mode_p) / cohort_spread[m]
            G[m] = v > gate
            levels.append({'image': t, 'marker': m, 'zero_fraction': float((x <= 0).mean()), 'log_mode': mode,
                           'log_spread': spread, 'background': float(np.expm1(mode_p)),
                           **{f'level_k{k}': float(np.expm1(mode + k * spread)) for k in KS},
                           'peak_mode': mode_p, 'peak_spread': spread_p,
                           **{f'peak_level_k{k}': float(np.expm1(mode_p + k * spread_p)) for k in KS},
                           'cohort_spread': cohort_spread[m],
                           **{f'cohort_level_k{k}': float(np.expm1(mode_p + k * cohort_spread[m])) for k in KS},
                           'scimap_gate': float(np.expm1(gate)),
                           'v15_level_linear': scale.loc[m, 'floor'] + 0.10 * scale.loc[m, 'span']})
        for k in KS:
            recs.append({'image': t, 'method': f'log, lower-half spread, z > {k}', **checks(Z > k)})
            recs.append({'image': t, 'method': f'log, peak-width spread, z > {k}', **checks(Zp > k)})
            recs.append({'image': t, 'method': f'log, cohort width, z > {k}', **checks(Zc > k)})
        recs.append({'image': t, 'method': 'scimap: log1p, GMM(2) gate', **checks(G)})
        # v15 for the same cells, on the aligned values the rules read
        al = pd.read_parquet(typ / 'aligned' / f'{t}.parquet', columns=['object_id'] + MARKERS)
        al = al.set_index('object_id').reindex(X['object_id'])
        S = pd.DataFrame(CL.apply_scale(al[MARKERS].to_numpy(), MARKERS, scale), columns=MARKERS)
        recs.append({'image': t, 'method': 'v15: linear scaled >= 0.10', **checks(S >= 0.10)})
        figure(t, X, {r['marker']: r for r in levels if r['image'] == t}, OUT / 'figures' / f'{t}_log_histograms.png')
        print(t, flush=True)
    R = pd.DataFrame(recs)
    R.to_csv(OUT / 'log_background_checks.csv', index=False)
    L = pd.DataFrame(levels)
    L.to_csv(OUT / 'log_background_levels.csv', index=False)
    pd.set_option('display.width', 300)
    keep = [m for m in R['method'].unique() if 'lower-half' not in m and 'z > 5' not in m]
    print(R[R['method'].isin(keep)].groupby('method').median(numeric_only=True).round(3).T.to_string())
    L['fold_peak_k3'] = L['peak_level_k3'] / L['background']
    L['fold_cohort_k3'] = L['cohort_level_k3'] / L['background']
    L['fold_v15'] = L['v15_level_linear'] / L['background']
    print('\ncohort width (log units) per marker:', {m: round(s, 3) for m, s in cohort_spread.items()})
    print('\nlevel as a multiple of the image background, median [min, max] over images:')
    for c in ['fold_v15', 'fold_peak_k3', 'fold_cohort_k3']:
        g = L.groupby('marker')[c].agg(['median', 'min', 'max']).round(1)
        print(c, {m: f'{r["median"]} [{r["min"]}, {r["max"]}]' for m, r in g.iterrows()})
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
