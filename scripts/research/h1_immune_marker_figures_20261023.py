#!/usr/bin/env python3
"""Scientific descriptive figures for all immune-associated and all31 markers.

Every heatmap cell is an equal-donor mean of a within-donor, unique-cell
summary. Gray represents unavailable measurement/support, not biological zero.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'islet-explorer-h1-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import h1_features_20261022 as hf
from h1_immune_marker_features_20261023 import MARKERS, IMMUNE_MARKERS, IMMUNE_TYPES, SIZE_BINS

STATUSES = ['ND', 'Aab+', 'T1D']
SIZE_LABELS = {'<50': '<50', '50-<100': '50–<100', '>=100': '≥100'}
TYPE_LABELS = {'CD8a Tcell': 'CD8 T', 'CD4 Tcell': 'CD4 T', 'T cell': 'T, untyped',
               'B cell': 'B', 'Macrophage': 'Macrophage', 'Immune': 'Other immune'}


def aggregate(d, group, value):
    eligible = d.loc[np.isfinite(d[value])].copy()
    return eligible.groupby(group, observed=True).agg(
        donor_balanced_mean=(value, 'mean'), donor_n=('case_id', 'nunique'),
        measured_cell_count=('n_measured', 'sum')).reset_index()


def heatmap(ax, data, rows, columns, label, fraction=True, show_y=True):
    matrix, ns = np.full((len(rows), len(columns)), np.nan), np.zeros((len(rows), len(columns)), int)
    lookup = data.set_index(['marker', 'column_key'])
    for i, marker in enumerate(rows):
        for j, key in enumerate(columns):
            if (marker, key) in lookup.index:
                r = lookup.loc[(marker, key)]
                matrix[i, j], ns[i, j] = r.donor_balanced_mean, int(r.donor_n)
    cmap = matplotlib.colormaps['YlGnBu'].copy()
    cmap.set_bad('#D5DADD')
    image = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=0, vmax=100 if fraction else 1, aspect='auto')
    for i in range(len(rows)):
        for j in range(len(columns)):
            v = matrix[i, j]
            number = '—' if not np.isfinite(v) else (f'{v:.0f}%' if fraction else f'{v:.3f}')
            color = 'white' if np.isfinite(v) and v > (57 if fraction else .57) else '#152F3D'
            ax.text(j, i, f'{number}\nn={ns[i, j]}', ha='center', va='center', fontsize=7 if fraction else 6.5,
                    linespacing=1.15, color=color)
    ax.set_yticks(range(len(rows)), rows if show_y else [''] * len(rows), fontsize=9 if fraction else 8.5)
    ax.set_xticks(range(len(columns)))
    ax.set_xticks(np.arange(-.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(rows), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.6)
    ax.tick_params(which='both', length=0)
    ax.set_title(label, fontsize=14, weight='bold', loc='left', pad=34 if fraction else 15, color='#173540')
    for spine in ax.spines.values():
        spine.set_visible(False)
    return image


def paired_size_figure(summary, target, title, subtitle, label):
    columns = [f'{status}|{size}' for status in STATUSES for size in SIZE_BINS]
    fig, axes = plt.subplots(1, 2, figsize=(18.5, 10))
    fig.subplots_adjust(left=.075, right=.895, top=.79, bottom=.17, wspace=.12)
    for k, compartment in enumerate(['core', 'peri']):
        d = summary.loc[summary.compartment.eq(compartment)].copy()
        d['column_key'] = d.donor_status + '|' + d.diameter_bin
        image = heatmap(axes[k], d, IMMUNE_MARKERS, columns, compartment.capitalize(), show_y=k == 0)
        axes[k].set_xticklabels([SIZE_LABELS[s] for _ in STATUSES for s in SIZE_BINS], fontsize=9)
        axes[k].set_xlabel('Equivalent core diameter (µm)', fontsize=10, labelpad=10)
        for g, status in enumerate(STATUSES):
            axes[k].text(g * 3 + 1, 1.025, status, transform=axes[k].get_xaxis_transform(), ha='center', fontsize=11, weight='bold')
        for v in [2.5, 5.5]:
            axes[k].axvline(v, color='#8D9DA6', lw=1.25)
    cb = fig.colorbar(image, cax=fig.add_axes([.92, .27, .012, .42]))
    cb.set_label(label, fontsize=10)
    fig.suptitle(title, x=.075, y=.965, ha='left', fontsize=21, weight='bold', color='#173540')
    fig.text(.075, .925, subtitle, fontsize=11, color='#526572')
    fig.text(.075, .083, 'Each cell is an equal-donor mean; n is donors with a measured, nonzero denominator. Gray: unavailable.',
             fontsize=10, color='#526572')
    fig.text(.075, .056, 'Cell IDs are deduplicated within donor × compartment × size bin; a shared peri cell can occur in different size bins.',
             fontsize=9.5, color='#526572')
    fig.text(.075, .029, 'Generic marker-positive profiles overlap. CD44, VIM, LGALS3 and other signals do not by themselves define an immune identity.',
             fontsize=9.5, color='#526572')
    for ext in ['png', 'svg']:
        fig.savefig(target.with_suffix('.' + ext), dpi=220)
    plt.close(fig)


def all31_figure(donor_table, output):
    d = donor_table.loc[donor_table.phenotype.isin(IMMUNE_TYPES) & donor_table.n_measured.gt(0)].copy()
    summary = aggregate(d, ['compartment', 'phenotype', 'marker'], 'mean')
    rows = IMMUNE_MARKERS + [m for m in MARKERS if m not in IMMUNE_MARKERS]
    fig, axes = plt.subplots(1, 2, figsize=(16, 17))
    fig.subplots_adjust(left=.09, right=.90, top=.90, bottom=.105, wspace=.13)
    for k, compartment in enumerate(['core', 'peri']):
        sub = summary.loc[summary.compartment.eq(compartment)].copy()
        sub['column_key'] = sub.phenotype
        image = heatmap(axes[k], sub, rows, IMMUNE_TYPES, compartment.capitalize(), fraction=False, show_y=k == 0)
        axes[k].set_xticklabels([TYPE_LABELS[x] for x in IMMUNE_TYPES], rotation=35, ha='right', fontsize=9)
        axes[k].axhline(11.5, color='#8D9DA6', lw=1.25)
    cb = fig.colorbar(image, cax=fig.add_axes([.925, .30, .012, .40]))
    cb.set_label('Mean scaled per-cell signal (0–1)', fontsize=10)
    fig.suptitle('All 31 markers within curated immune types', x=.09, y=.975, ha='left', fontsize=21, weight='bold', color='#173540')
    fig.text(.09, .947, 'Equal-donor arithmetic means after deduplicating local cells within donor and compartment; no marker-based retyping.',
             fontsize=10.5, color='#526572')
    fig.text(.09, .055, 'n = donors with ≥1 measured cell of the indicated type; sparse core types remain descriptive. Gray = no measured donor support.',
             fontsize=9.5, color='#526572')
    fig.text(.09, .034, 'Cytoplasmic measurements are unavailable in some nucleus-fills-cell cases; means use measured cells only. No missing-value imputation.',
             fontsize=9.5, color='#526572')
    fig.text(.09, .014, 'Scaled fluorescence is not calibrated protein concentration. Typing-marker patterns partly reflect the gates that named each cell.',
             fontsize=9.5, color='#526572')
    for ext in ['png', 'svg']:
        fig.savefig(output.with_suffix('.' + ext), dpi=220)
    plt.close(fig)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feature-root', type=Path, required=True)
    args = parser.parse_args()
    root = args.feature_root.resolve(strict=True)
    out = root / 'figures'
    out.mkdir(exist_ok=True)
    sizes = pd.read_csv(root / 'donor_size_type_compartment_markers.csv')
    donors = pd.read_csv(root / 'donor_type_compartment_markers.csv')
    immune = sizes.loc[sizes.phenotype.eq('IMMUNE_UNION') & sizes.marker.isin(IMMUNE_MARKERS)].copy()
    immune['percentage'] = np.where(immune.n_measured.gt(0), 100 * immune.n_positive / immune.n_measured, np.nan)
    summary = aggregate(immune, ['compartment', 'donor_status', 'diameter_bin', 'marker'], 'percentage')
    summary.to_csv(out / 'fig3_immune_marker_profiles_data.csv', index=False)
    paired_size_figure(summary, out / 'fig3_immune_marker_profiles', 'Immune-associated marker profiles by size and donor status',
                       'Marker-positive percentage within the curated immune union • inclusive scaled signal ≥0.10 • all 12 requested markers',
                       'Positive among measured immune cells (%)')
    all_cells = sizes.loc[sizes.phenotype.eq('ALL') & sizes.marker.isin(IMMUNE_MARKERS)].copy()
    keys = ['case_id', 'donor_status', 'compartment', 'diameter_bin', 'marker']
    shared = all_cells.merge(immune[keys + ['n_cells', 'n_measured', 'n_positive']], on=keys,
                             suffixes=('_all', '_immune'), validate='one_to_one')
    available = shared.n_positive_all.gt(0) & shared.n_measured_all.gt(0) & (
        shared.n_cells_immune.eq(0) | shared.n_measured_immune.gt(0))
    shared['percentage_outside_immune'] = np.where(available,
        100 * (shared.n_positive_all - shared.n_positive_immune) / shared.n_positive_all, np.nan)
    shared['n_measured'] = shared.n_positive_all
    outside_summary = aggregate(shared, ['compartment', 'donor_status', 'diameter_bin', 'marker'], 'percentage_outside_immune')
    outside_summary.to_csv(out / 'figS2_marker_positive_outside_immune_data.csv', index=False)
    paired_size_figure(outside_summary, out / 'figS2_marker_positive_outside_immune', 'Marker positivity outside the curated immune union',
                       'Among all measured marker-positive local cells, percentage lacking a curated immune identity • inclusive scaled signal ≥0.10',
                       'Marker-positive cells outside immune union (%)')
    all_summary = all31_figure(donors, out / 'figS1_all31_immune_profiles')
    all_summary.to_csv(out / 'figS1_all31_immune_profiles_data.csv', index=False)
    audit = dict(aggregation='equal donor, unique cells within donor+compartment or donor+compartment+size',
                 no_imputation=True, primary_marker_threshold=.10, primary_markers=IMMUNE_MARKERS,
                 all31_markers=MARKERS, statuses=STATUSES, diameter_bins=SIZE_BINS,
                 denominator='measured immune cells for fig3; observed marker-positive all cells for figS2; measured type cells for figS1',
                 missing_color='#D5DADD', source_sha256={str(root / n): hf.sha256(root / n) for n in
                    ['donor_size_type_compartment_markers.csv', 'donor_type_compartment_markers.csv']},
                 script_sha256=hf.sha256(Path(__file__)),
                 output_sha256={str(p): hf.sha256(p) for p in out.glob('*.png') if p.name.startswith(('fig3_', 'figS1_', 'figS2_'))},
                 descriptive_observations={
                    'median_outside_immune_percentage_by_marker': outside_summary.groupby('marker').donor_balanced_mean.median().to_dict(),
                    'fig3_donor_support_range': [int(summary.donor_n.min()), int(summary.donor_n.max())],
                    'all31_core_type_donor_support_range': [int(all_summary.loc[all_summary.compartment.eq('core')].donor_n.min()),
                                                          int(all_summary.loc[all_summary.compartment.eq('core')].donor_n.max())]})
    hf.write_json(out / 'marker_figure_audit.json', audit)
    print(json.dumps(hf.finite_json(audit['descriptive_observations']), indent=2))


if __name__ == '__main__':
    main()
