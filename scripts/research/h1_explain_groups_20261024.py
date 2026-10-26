#!/usr/bin/env python3
"""Explain saved corrected H1 estimates; no models are refitted."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/islet-h1-matplotlib')
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'output/h1_edge_corrected_20261024'
OUT = ROOT / 'output/h1_methods_review_20261024'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    models = pd.read_csv(SOURCE / 'models/model_results.csv').set_index('model')
    slopes = pd.read_csv(SOURCE / 'heterogeneity/group_specific_slopes.csv')
    d = pd.read_csv(SOURCE / 'models/primary_analysis_objects.csv')
    records = []
    for donor, sub in d[d.donor_status.eq('T1D')].groupby('case_id'):
        for inside, cells in sub.groupby('has_core_immune'):
            records.append(dict(case_id=int(donor), core_immune_present=int(inside),
                objects=len(cells), beta_present=int(cells.n_beta_core.gt(0).sum()),
                percent_beta_present=100*cells.n_beta_core.gt(0).mean()))
    census = pd.DataFrame(records)
    census.to_csv(OUT / 't1d_beta_presence.csv', index=False)
    assert census.objects.sum() == 1088 and census.beta_present.sum() == 80
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':11,
                         'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios':[1.1,1]})
    rows = []
    for name, label in [('Primary','Pooled: common slope'), ('T1D only','T1D: separate model')]:
        r = models.loc[name]
        rows.append((label,r.rr,r.lower,r.upper,'#26364A'))
    for r in slopes.itertuples():
        rows.append((f'{r.donor_status}: group slope*',r.rr_10pp,r.lower_95,r.upper_95,
                     {'ND':'#4477AA','Aab+':'#CC6633','T1D':'#228833'}[r.donor_status]))
    for i,(label,rr,lo,hi,color) in enumerate(rows):
        axes[0].errorbar(rr,i,xerr=[[rr-lo],[hi-rr]],fmt='o',color=color,capsize=4)
    axes[0].axvline(1,color='gray',ls='--',lw=1)
    axes[0].axhline(1.5,color='#ddd',lw=1)
    axes[0].set_yticks(range(len(rows)),[r[0] for r in rows])
    axes[0].invert_yaxis()
    axes[0].set_xlabel('Beta-density ratio per 10 percentage points\n95% intervals; values below 1 are inverse')
    axes[0].set_title('A  A common slope hides different group patterns',loc='left',fontsize=12)
    colors = plt.cm.tab10(np.linspace(0, .7, 6))
    for color,(donor,sub) in zip(colors,census.groupby('case_id')):
        sub = sub.sort_values('core_immune_present')
        axes[1].plot([0,1],sub.percent_beta_present,marker='o',color=color,lw=1.6,label=str(donor))
    axes[1].set_xticks([0,1],['Immune cells\naround only','Immune cells\nalso in core'])
    axes[1].set_xlim(-.12,1.12)
    axes[1].set_ylim(0,70)
    axes[1].set_ylabel('Objects with ≥1 typed Beta cell (%)')
    axes[1].set_title('B  Residual Beta cells and core occupancy in T1D',loc='left',fontsize=12)
    axes[1].legend(title='Donor',ncol=2,fontsize=9,loc='upper left')
    fig.text(.02,.025,'*Posthoc group slopes share nuisance coefficients across all donors; approximate t(19) intervals.\n'
             'Panel B is unadjusted and descriptive: all six donors rise, but object size and other covariates are not controlled.',fontsize=9)
    fig.tight_layout(rect=[0,.12,1,1])
    fig.savefig(OUT/'pooled_vs_t1d.png',dpi=180)
    fig.savefig(OUT/'pooled_vs_t1d.pdf')
    plt.close(fig)


if __name__ == '__main__':
    main()
