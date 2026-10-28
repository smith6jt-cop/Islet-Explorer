#!/usr/bin/env python3
"""Is the immune positivity level (scaled >= 0.10) too high? Owner request, 2026-10-27: "what is clear now is that the
top marker of leftover cells indicates the threshold is too high for positivity in many cases. Investigate and report."

For every marker the typing reads, on root v15 and the run's aligned values and cohort scale:
  - the 0.10 level in raw (aligned) units;
  - its background in cells that should lack it (negative reference): endocrine cells (Beta / Alpha / Delta, typed
    by hormone inside the endocrine objects) and rule-typed Exocrine cells for immune, stromal and vascular markers;
    CD3e-high cells (scaled CD3e >= 0.5, selected on CD3e alone) for E-cadherin and panCK; immune rule-typed cells for
    the hormones. Reported: the reference's p95 / p99 on the cohort scale (before the clip) and the proportion at or
    above 0.10;
  - an independent positive reference where one exists: CD45 in CD3e-high cells;
  - the left-over cells the top-marker step placed or named by this marker: their value on it, and the proportion
    above the reference's p99 and above 0.10.
Per image, the same for the immune markers. Measurements only.

    python scripts/research/positivity_threshold_review_20261027.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
ROOT = REPO / 'data' / 'development' / 'uniform-20donor-20261027-v15'
OUT = REPO / 'output' / 'positivity_threshold_review_20261027'
LEVEL = CT.IMMUNE_REQUIRE_MIN
ENDOCRINE = ('Beta cell', 'Alpha cell', 'Delta cell')
IMMUNE_RULE = ('CD8a Tcell', 'CD4 Tcell', 'T cell', 'B cell', 'Immune', 'CD45+')
# marker -> which negative reference applies
NEGATIVE = {**{m: 'endocrine + rule Exocrine' for m in ('CD3e', 'CD8a', 'CD4', 'CD20', 'CD45', 'CD68', 'CD163',
                                                         'LGALS3', 'CD31', 'CD34', 'PDPN', 'VIM', 'SMA', 'B3TUBB',
                                                         'ColIV')},
            'ECAD': 'CD3e-high cells', 'panCK': 'CD3e-high cells',
            'INS': 'immune rule-typed', 'GCG': 'immune rule-typed', 'SST': 'immune rule-typed'}


def load() -> pd.DataFrame:
    typ = CT.typing_dir(RUN)
    frames = []
    markers = sorted(NEGATIVE)
    for r in CT.cohort_rows():
        t = r['image_name']
        lab = pd.read_parquet(ROOT / 'cell_types' / f"{int(r['case_id'])}.parquet", columns=['object_id', 'phenotype'])
        calls = pd.read_parquet(typ / 'calls' / f'{t}.parquet', columns=['object_id', 'resolved_by'])
        al = pd.read_parquet(typ / 'aligned' / f'{t}.parquet', columns=['object_id'] + markers)
        d = lab.merge(calls, on='object_id').merge(al, on='object_id')
        d['image'] = t
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    scale = CL.load_scale(RUN)
    d = load()
    markers = sorted(NEGATIVE)
    U = pd.DataFrame({m: (d[m] - scale.loc[m, 'floor']) / scale.loc[m, 'span'] for m in markers})   # before the clip
    lab = d['phenotype']
    top_rows = d['resolved_by'].str.endswith('top_marker').to_numpy()
    # which marker each left-over cell was placed or named by: its top marker over the profile, as classify did
    prof = [m for m in CT.TOP_MARKER_PROFILE if m in d.columns]
    Up = pd.DataFrame({m: (d[m] - scale.loc[m, 'floor']) / scale.loc[m, 'span'] for m in prof})
    top_marker = Up.idxmax(axis=1).where(top_rows, '')
    cd3e_high = U['CD3e'] >= 0.5
    refs = {'endocrine + rule Exocrine': lab.isin(ENDOCRINE) | ((lab == 'Exocrine') & (d['resolved_by'] == 'definition_exocrine')),
            'CD3e-high cells': cd3e_high, 'immune rule-typed': lab.isin(IMMUNE_RULE) | lab.str.contains('Macrophage_')}
    rows = []
    for m in markers:
        ref = refs[NEGATIVE[m]].to_numpy()
        neg = U.loc[ref, m].to_numpy()
        neg = neg[np.isfinite(neg)]
        p95, p99 = (np.percentile(neg, [95, 99]) if len(neg) else (np.nan, np.nan))
        lo = U.loc[(top_marker == m).to_numpy(), m].to_numpy()
        rows.append({'marker': m, 'floor_raw': scale.loc[m, 'floor'], 'top_raw': scale.loc[m, 'top'],
                     'level_0.10_raw': scale.loc[m, 'floor'] + LEVEL * scale.loc[m, 'span'],
                     'negative_reference': NEGATIVE[m], 'negative_n': len(neg),
                     'negative_p95_scaled': p95, 'negative_p99_scaled': p99,
                     'negative_p99_raw': scale.loc[m, 'floor'] + p99 * scale.loc[m, 'span'],
                     'negative_proportion_at_or_above_0.10': float((neg >= LEVEL).mean()) if len(neg) else np.nan,
                     'level_over_negative_p99': LEVEL / p99 if p99 and p99 > 0 else np.nan,
                     'leftover_by_this_top_marker': len(lo),
                     'leftover_median_scaled': float(np.median(lo)) if len(lo) else np.nan,
                     'leftover_proportion_above_negative_p99': float((lo > p99).mean()) if len(lo) else np.nan,
                     'leftover_proportion_at_or_above_0.10': float((lo >= LEVEL).mean()) if len(lo) else np.nan})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / 'positivity_by_marker.csv', index=False)

    # independent positive reference: CD45 in CD3e-high cells
    pos = U.loc[cd3e_high, 'CD45'].to_numpy()
    pos_rec = {'marker': 'CD45', 'positive_reference': 'CD3e-high cells (scaled CD3e >= 0.5)', 'n': int(len(pos)),
               'median_scaled': float(np.median(pos)), 'proportion_at_or_above_0.10': float((pos >= LEVEL).mean()),
               'proportion_above_negative_p99': float((pos > tab.set_index('marker').loc['CD45', 'negative_p99_scaled']).mean())}
    pd.DataFrame([pos_rec]).to_csv(OUT / 'positive_reference_CD45.csv', index=False)

    # per image, the immune markers
    per = []
    for t, g in d.groupby('image'):
        Ug = U.loc[g.index]
        for m in ('CD3e', 'CD8a', 'CD4', 'CD20', 'CD45', 'CD68', 'CD163'):
            neg = Ug.loc[refs[NEGATIVE[m]].loc[g.index].to_numpy(), m].to_numpy()
            neg = neg[np.isfinite(neg)]
            per.append({'image': t, 'marker': m, 'negative_n': len(neg),
                        'negative_p99_scaled': float(np.percentile(neg, 99)) if len(neg) else np.nan,
                        'negative_proportion_at_or_above_0.10': float((neg >= LEVEL).mean()) if len(neg) else np.nan})
    pd.DataFrame(per).to_csv(OUT / 'negative_p99_by_image.csv', index=False)

    # figures: negative reference vs left-over cells on their top marker, scaled before the clip
    show = [m for m in ('CD20', 'CD45', 'CD4', 'CD8a', 'CD3e', 'CD68', 'CD163', 'LGALS3', 'CD34', 'VIM', 'CD31', 'ECAD')
            if m in markers]
    fig, axes = plt.subplots(3, 4, figsize=(20, 13))
    for ax, m in zip(axes.flat, show):
        neg = U.loc[refs[NEGATIVE[m]].to_numpy(), m].dropna().to_numpy()
        lo = U.loc[(top_marker == m).to_numpy(), m].dropna().to_numpy()
        bins = np.linspace(-0.05, 0.4, 91)
        ax.hist(np.clip(neg, bins[0], bins[-1]), bins=bins, density=True, alpha=0.55, color='#4477AA',
                label=f'negative reference (n={len(neg):,})')
        if len(lo):
            ax.hist(np.clip(lo, bins[0], bins[-1]), bins=bins, density=True, alpha=0.55, color='#CC6633',
                    label=f'left-over cells, top marker {m} (n={len(lo):,})')
        p99 = tab.set_index('marker').loc[m, 'negative_p99_scaled']
        ax.axvline(LEVEL, color='black', lw=1.5, label='0.10 level')
        ax.axvline(p99, color='#228833', lw=1.5, ls='--', label='negative p99')
        ax.set_title(m, fontsize=15)
        ax.set_xlabel('cohort-scaled value (before the 0-1 clip)', fontsize=12)
        ax.legend(fontsize=9)
    fig.suptitle('Marker values: cells that should lack the marker vs left-over cells placed or named by it',
                 fontsize=17)
    fig.tight_layout()
    fig.savefig(OUT / 'negative_vs_leftover_histograms.png', dpi=150)
    print(tab[['marker', 'level_0.10_raw', 'negative_p99_scaled', 'negative_proportion_at_or_above_0.10',
               'leftover_by_this_top_marker', 'leftover_median_scaled', 'leftover_proportion_above_negative_p99',
               'leftover_proportion_at_or_above_0.10']].round(3).to_string(index=False))
    print(pos_rec)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
