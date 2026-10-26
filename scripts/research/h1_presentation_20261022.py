#!/usr/bin/env python3
"""Plot stored H1 estimates without selecting or refitting statistical models."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/islet-h1-matplotlib')
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path('output/h1_20261022')
GROUPS = ['ND', 'Aab+', 'T1D']
COLORS = {'ND': '#4477AA', 'Aab+': '#CC6633', 'T1D': '#228833'}


def save(fig, name):
    fig.savefig(ROOT / (name+'.png'), dpi=220, bbox_inches='tight', facecolor='white')
    fig.savefig(ROOT / (name+'.svg'), bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titlesize': 12, 'axes.labelsize': 11})
    groups = pd.read_csv(ROOT/'heterogeneity/group_specific_slopes.csv')
    models = pd.read_csv(ROOT/'model_results.csv')
    matched = pd.read_csv(ROOT/'independent/donor_effects.csv', dtype={'case_id': str})
    matched = matched.loc[matched.n_pairs >= 5].copy()
    matched['order'] = matched.donor_status.map(dict(zip(GROUPS, range(3))))
    matched = matched.sort_values(['order', 'case_id'])
    ms = json.loads((ROOT/'independent/summary.json').read_text())
    primary = models.loc[models.model == 'Primary'].iloc[0]
    t1d = models.loc[models.model == 'T1D only'].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(13, 7.5), gridspec_kw={'width_ratios': [1.05, 1.25]})
    axes[0].axvline(1, color='#999999', ls='--', lw=1)
    ys = [0, 1.5, 2.5, 3.5, 5]
    labels = ['All donors: primary', 'ND: shared-control model',
              'Aab+: shared-control model', 'T1D: shared-control model',
              'T1D: separate sensitivity']
    rr = [primary.rr]+groups.rr_10pp.tolist()+[t1d.rr]
    lo = [primary.lower]+groups.lower_95.tolist()+[t1d.lower]
    hi = [primary.upper]+groups.upper_95.tolist()+[t1d.upper]
    cs = ['#26364A']+[COLORS[g] for g in GROUPS]+[COLORS['T1D']]
    for y, v, a, b, c in zip(ys, rr, lo, hi, cs):
        axes[0].errorbar(v, y, xerr=[[v-a], [b-v]], fmt='o', color=c, capsize=4, markersize=7)
    axes[0].set_yticks(ys, labels)
    axes[0].invert_yaxis()
    axes[0].set_xscale('log')
    ticks = [.8, .9, 1, 1.1, 1.2, 1.4]
    axes[0].xaxis.set_major_locator(FixedLocator(ticks))
    axes[0].xaxis.set_major_formatter(FixedFormatter([f'{x:.2f}' for x in ticks]))
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].set_xlabel('Beta-density ratio per 10 percentage points inside\nBelow 1: lower density; above 1: higher density')
    axes[0].set_title('A  The direction differs across donor groups', loc='left', fontweight='bold')
    y = np.arange(len(matched))
    for j, row in enumerate(matched.itertuples()):
        axes[1].scatter(row.beta_density_difference_mean, j, color=COLORS[row.donor_status], s=45, zorder=3)
    axes[1].set_yticks(list(y)+[len(y)+1],
        [f'{r.case_id} {r.donor_status} ({r.n_pairs} pairs)' for r in matched.itertuples()]+['Equal-donor mean (12 donors)'])
    mean = ms['donor_balanced_beta_density_difference_cells_per_1000um2']
    a, b = ms['bootstrap_95pct_interval']
    axes[1].errorbar(mean, len(y)+1, xerr=[[mean-a], [b-mean]], fmt='D', color='#26364A', capsize=4, markersize=8)
    axes[1].axvline(0, color='#999999', ls='--', lw=1)
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Matched Beta-density difference (cells per 1,000 µm²)\nCore immune present minus absent')
    axes[1].set_title('B  Matched objects also show donor heterogeneity', loc='left', fontweight='bold')
    fig.suptitle('A pooled inverse association does not describe every donor or group', fontsize=15, fontweight='bold')
    fig.text(.5, .015, 'A: donor-cluster 95% intervals. Group-specific shared-control slopes are posthoc; the separate T1D-only sensitivity was specified before fitting.\n'
             'Only six T1D donors: subgroup estimates are exploratory. B: points are donor means; the pooled diamond has a donor-bootstrap interval.\n'
             'Matching uses the same dataset and a different, binary exposure; it is a sensitivity check, not independent replication.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .11, 1, .94])
    save(fig, 'figure4_groups_and_matching')

    # All eligible objects remain in the census; no Beta-based selection.
    census = pd.read_csv(ROOT/'heterogeneity/four_state_census_by_donor.csv', dtype={'case_id': str})
    support = pd.read_csv(ROOT/'donor_analysis_support.csv', dtype={'case_id': str})
    states = [(True, False), (True, True), (False, True), (False, False)]
    palette = ['#E7A58E', '#BF5634', '#377EAC', '#D7DCE2']
    names = ['Beta present / no core immune', 'Beta present / core immune',
             'No typed Beta / core immune', 'No typed Beta / no core immune']
    fig, axes = plt.subplots(1, 2, figsize=(13, 8.5), gridspec_kw={'width_ratios': [1.4, 1]})
    ys = np.arange(len(support))
    left = np.zeros(len(support))
    for state, color, name in zip(states, palette, names):
        s = census.loc[(census.beta_present == state[0]) & (census.core_immune_present == state[1])].set_index('case_id')
        vals = s.reindex(support.case_id).proportion_among_exposed.to_numpy()*100
        axes[0].barh(ys, vals, left=left, color=color, height=.75, label=name)
        left += vals
    axes[0].set_yticks(ys, [f'{r.case_id} {r.donor_status}' for r in support.itertuples()])
    axes[0].invert_yaxis()
    axes[0].set_xlim(0,100)
    axes[0].set_xlabel('Eligible nearby-immune-exposed objects (%)')
    axes[0].set_title('A  Current cell states vary within each donor', loc='left', fontweight='bold')
    axes[0].legend(loc='lower center', bbox_to_anchor=(.5,-.22), fontsize=9, ncol=2)
    two = models.loc[models.model == 'Two-part localization'].copy()
    for j, r in enumerate(two.itertuples()):
        axes[1].errorbar(r.rr, j, xerr=[[r.rr-r.lower], [r.upper-r.rr]], fmt='o', color='#26364A', capsize=4)
    axes[1].axvline(1, color='gray', ls='--')
    axes[1].set_yticks([0,1], ['Any core immune cells\nat median positive allocation',
                              'Additional 10 percentage points\namong core-positive objects'])
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Adjusted Beta-density ratio (95% CI)')
    axes[1].set_title('B  Separating presence from allocation dose', loc='left', fontweight='bold')
    axes[1].set_ylim(1.7,-.7)
    axes[1].text(.04,.08,'Both intervals include 1.\nThe data do not establish a distinct\nBeta-density drop at first core occupancy.',
                  transform=axes[1].transAxes, va='bottom', fontsize=10,
                  bbox={'boxstyle':'round,pad=.5','facecolor':'#F2F5F8','edgecolor':'#D5DDE5'})
    fig.suptitle('A snapshot of local immune occupancy is not a record of prior attack', fontsize=15, fontweight='bold')
    fig.text(.5,.015,'Beta absence means no typed Beta cells in that analysed core; it does not establish cell death, origin, or past immune exposure.\n'
             'The census is descriptive. Two-part estimates are the prespecified secondary model, with donor-cluster uncertainty.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0,.09,1,.94])
    save(fig, 'figure6_states_and_two_part')

    # Posthoc target-availability description: show the six independent donors.
    t = census.loc[census.donor_status == 'T1D'].copy()
    rows = []
    for case, sub in t.groupby('case_id'):
        row = {'case_id': case}
        for present, key in [(False, 'around_only'), (True, 'core_present')]:
            ss = sub.loc[sub.core_immune_present == present]
            num = int(ss.loc[ss.beta_present, 'objects'].sum())
            den = int(ss.objects.sum())
            row[key+'_beta_present'] = num
            row[key+'_objects'] = den
            row[key+'_percent'] = 100*num/den
        rows.append(row)
    td = pd.DataFrame(rows)
    td.to_csv(ROOT/'t1d_beta_presence_by_core_occupancy.csv', index=False)
    fig, ax = plt.subplots(figsize=(10,5.5))
    for j, r in enumerate(td.itertuples()):
        x0, x1 = r.around_only_percent, r.core_present_percent
        ax.plot([x0,x1], [j,j], color='#AEB9B1', lw=2, zorder=1)
        ax.scatter(x0,j,color='#7C8792',s=65,zorder=3)
        ax.scatter(x1,j,color=COLORS['T1D'],s=75,zorder=3)
        ax.annotate(f'{r.around_only_beta_present}/{r.around_only_objects}', (x0,j),
                    xytext=(0,10), textcoords='offset points', ha='center', fontsize=9, color='#5B6570')
        ax.annotate(f'{r.core_present_beta_present}/{r.core_present_objects}', (x1,j),
                    xytext=(0,-17), textcoords='offset points', ha='center', fontsize=9, color='#17633A')
    ax.set_yticks(range(len(td)), ['T1D '+str(x) for x in td.case_id])
    ax.set_ylim(len(td)-.35,-.55)
    ax.set_xlim(0,65)
    ax.set_xlabel('Objects containing at least one typed Beta cell (%)')
    ax.set_title('In all six T1D donors, objects with core immune cells\n'
                 'more often contain typed Beta cells', loc='left', fontsize=14, fontweight='bold', pad=22)
    ax.scatter([],[],color='#7C8792',s=65,label='Immune cells around only')
    ax.scatter([],[],color=COLORS['T1D'],s=75,label='Some immune cells inside core')
    ax.legend(loc='lower right',fontsize=10)
    fig.text(.5,.015,'Labels: Beta-containing objects / objects in that occupancy category. Unadjusted posthoc description; object size is not controlled here.\n'
             'Pooled counts: 51/168 with core immune cells versus 55/1,058 with immune cells around only. No object-level significance test.\n'
             'Size, composition and target availability are competing explanations; this does not establish antigen specificity or a temporal sequence.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.13,1,1])
    save(fig,'figure7_t1d_beta_presence')


if __name__ == '__main__':
    main()
