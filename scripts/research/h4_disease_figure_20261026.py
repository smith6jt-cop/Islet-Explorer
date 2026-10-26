#!/usr/bin/env python3
"""Disease-state standardized regression coefficients and 95% confidence intervals."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h4-disease-20261026-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT=Path('output/h1_h4_disease_comparisons_20261026/h4')
GROUPS=['ND','Aab+','T1D'];COLORS=['#4477AA','#CC6633','#228833']

def main():
    d=pd.read_csv(OUT/'disease_slopes.csv')
    d=d.loc[d.model.eq('primary')].set_index('donor_status').loc[GROUPS]
    p=pd.read_csv(OUT/'disease_slope_comparisons.csv')
    p=p.loc[p.model.eq('primary')]
    a=pd.read_csv(OUT/'model_anova.csv')
    a=a.loc[a.model.eq('primary')&a.term.eq('donor_status:source_z')].iloc[0]
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(6.7,5.1))
    ax.axhline(0,lw=1,color='#999999',linestyle='--',zorder=1)
    for i,(_,r) in enumerate(d.iterrows()):
        ax.errorbar(i,r['source_z.trend'],
            yerr=[[r['source_z.trend']-r['lower.CL']],[r['upper.CL']-r['source_z.trend']]],
            fmt='o',color=COLORS[i],ms=9,lw=2,capsize=6,zorder=3)
    ax.set_xticks(range(3),GROUPS)
    ax.set_xlim(-.55,2.55)
    ax.set_ylim(min(d['lower.CL'])-.025,max(d['upper.CL'])+.025)
    ax.set_ylabel('Standardized regression coefficient (β)')
    ax.set_title('LGALS3 source allocation and β-cell density',fontsize=13,pad=17)
    # Comparison brackets are only shown for BH-significant disease contrasts.
    # None are significant in the primary fitted analysis.
    assert not p.p_bh.lt(.05).any()
    fig.tight_layout()
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'source_association_by_disease.{ext}',dpi=200)
    plt.close(fig)
    rows=[]
    for group,r in d.iterrows():
        rows.append(f'| {group} | {int(r["donors"])} | {int(r["objects"]):,} | '
                    f'{r["source_z.trend"]:.3f} [{r["lower.CL"]:.3f}, {r["upper.CL"]:.3f}] |')
    txt=(
      '# Figure legend\n\n'
      'LGALS3 source association by disease state. Points show adjusted standardized '
      'regression coefficients for myeloid LGALS3 source share versus Beta-cell density; '
      'bars show 95% confidence intervals. These are regression coefficients, not '
      'correlation coefficients. Source share and continuous covariates were standardized '
      'within donor, and Beta density was standardized across the included endocrine objects. '
      'The model contains donor random intercepts and source slopes, disease-by-source '
      'interaction, total source signal, myeloid and epithelial cell densities, immune '
      'density, core area, and object class.\n\n'
      f'Disease-by-source interaction: F(2, {a["DenDF"]:.2f})={a["F value"]:.3f}, '
      f'p={a["Pr(>F)"]:.3f}. BH-adjusted pairwise disease-comparison p values are '
      + '; '.join(f'{r["contrast"]}: {r["p_bh"]:.3f}' for _,r in p.iterrows())
      + '. No pairwise comparison is significant, so no significance brackets are shown.\n\n'
      '| Disease state | Donors | Endocrine objects | Standardized β [95% CI] |\n'
      '|---|---:|---:|---:|\n'+'\n'.join(rows)+'\n')
    (OUT/'figure_legend.md').write_text(txt)
    print(txt)

if __name__=='__main__':main()
