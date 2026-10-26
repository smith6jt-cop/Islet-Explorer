#!/usr/bin/env python3
"""Render the saved donor-level disease comparisons without refitting models."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h1-disease-figure-20261026')
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h1_h4_disease_comparisons_20261026/h1'
GROUPS=['ND','Aab+','T1D']
COLORS=['#4477AA','#CC6633','#228833']


def main():
    d=pd.read_csv(OUT/'donor_coefficients.csv')
    comparisons=pd.read_csv(OUT/'donor_pairwise_disease_tests.csv')
    overall=pd.read_csv(OUT/'donor_disease_global_test.csv').iloc[0]
    samples=[d.loc[d.donor_status.eq(g),'standardized_beta'].to_numpy() for g in GROUPS]
    fig,ax=plt.subplots(figsize=(6.2,5.1))
    boxes=ax.boxplot(samples,positions=np.arange(3),widths=.48,patch_artist=True,
                     showfliers=False,medianprops={'color':'#222222','linewidth':1.5})
    for i,(values,color,box) in enumerate(zip(samples,COLORS,boxes['boxes'])):
        box.set_facecolor(color);box.set_alpha(.18);box.set_edgecolor(color)
        ax.scatter(i+np.linspace(-.09,.09,len(values)),values,s=45,color=color,
                   edgecolor='white',linewidth=.6,zorder=3)
    ax.axhline(0,color='#888888',ls='--',lw=.8,zorder=0)
    ax.set_xticks(np.arange(3),GROUPS)
    ax.set_xlabel('Disease state')
    ax.set_ylabel('Adjusted association with β-cell density\n(standardized regression coefficient)')
    ax.set_title('Immune localization and β-cell density',fontsize=13,pad=15)
    ax.spines[['top','right']].set_visible(False)
    low,high=float(d.standardized_beta.min()),float(d.standardized_beta.max())
    span=high-low
    top=high
    for r in comparisons.loc[comparisons.p_bh.lt(.05)].itertuples():
        i,j=GROUPS.index(r.group_1),GROUPS.index(r.group_2)
        base=top+.10*span;height=.03*span
        ax.plot([i,i,j,j],[base,base+height,base+height,base],color='black',lw=1)
        star='***' if r.p_bh<.001 else '**' if r.p_bh<.01 else '*'
        ax.text((i+j)/2,base+height,star,ha='center',va='bottom',fontsize=16)
        top+=.20*span
    ax.set_ylim(low-.10*span,top+.12*span)
    fig.tight_layout()
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/f'immune_beta_association_by_disease.{ext}',dpi=220)
    plt.close(fig)
    pair_text='; '.join(f'{r.group_1} versus {r.group_2}, adjusted p={r.p_bh:.4g}'
                       for r in comparisons.itertuples())
    legend=(
        'Immune localization and β-cell density across disease states. Each point is one donor’s '
        'covariate-adjusted regression coefficient relating the fraction of nearby immune cells inside '
        'the endocrine core to core β-cell density. The nearby immune pool is the unique union '
        'of cells in the core and its surrounding 20-µm peri band. Positive values indicate higher β-cell density '
        'with greater core immune allocation; negative values indicate lower density. These are '
        'regression coefficients, not correlation coefficients. The regressions adjust for core area, '
        'peri-to-core area ratio, nearby immune-cell count, immune density in the independent '
        '20–100-µm surrounding field, and '
        'islet/cluster class. Exposure and continuous covariates were standardized within donor; '
        'β-cell density was standardized across all 2,937 eligible endocrine objects. Thus the '
        'coefficient expresses the adjusted change in outcome standard deviations for one '
        'within-donor exposure standard deviation. ND, n=7 donors (792 objects); Aab+, n=7 '
        '(1,057 objects); T1D, n=6 (1,088 objects). Boxes show medians and interquartile ranges; '
        'whiskers extend to the most extreme values within 1.5 interquartile ranges. All donor '
        'values are shown. Disease groups were compared using a Kruskal–Wallis test '
        f'(H={overall.statistic:.3f}, df=2, p={overall.p_value:.4g}), followed by all three '
        'two-sided exact Wilcoxon rank-sum tests with Benjamini–Hochberg adjustment: '
        f'{pair_text}. The bracket marks the significant ND–T1D pair (*adjusted p<0.05). '
        'This adjustment is across the three pairwise comparisons within this hypothesis. '
        'The separate adjustment across all four overall hypothesis tests is reported in the '
        'combined results. ND, non-diabetic; Aab+, autoantibody-positive; T1D, type 1 diabetes.'
    )
    (OUT/'figure_legend.txt').write_text(legend+'\n')
    report=(
        '# Immune localization–β-cell relationship across disease states\n\n'
        f'The direct donor-level comparison gives Kruskal–Wallis p={overall.p_value:.5f}. '
        'The ND–T1D comparison has BH-adjusted p=.04196; the ND–Aab+ and Aab+–T1D '
        'comparisons have adjusted p=.10927 and .18065. Median standardized adjusted '
        'coefficients are −.1274 (ND), −.0116 (Aab+), and +.0176 (T1D).\n\n'
        '![Disease comparison](immune_beta_association_by_disease.png)\n\n'+legend+'\n\n'
        'The requested donor random-intercept/random-slope models were fitted first. Both '
        'converged numerically but had singular random-effect covariance: the correlated '
        'model reached an intercept–slope correlation of −1, and the independent-slope '
        'model estimated zero slope variance. Their tests are preserved as diagnostics. '
        'The primary donor comparison therefore uses one separately adjusted coefficient '
        'per donor, not the object-scale degrees of freedom from those boundary fits. '
        'All 20 separate regressions have full rank.\n\n'
        'Methods: Barlow PDF p22 supports within-donor predictor standardization, global '
        'outcome standardization and hierarchical linear models. Richardson PDF p8 supports '
        'Kruskal–Wallis disease comparisons; Bodenmiller PDF p28 supports donor-aggregated '
        'Wilcoxon rank-sum comparisons and FDR adjustment. Small untied samples use the exact '
        'rank-sum calculation; tied pairs use a tie-corrected approximation with continuity '
        'correction. All three H1 pairs are untied.\n\n'
        'Source data and the previous H1 results are preserved. Independent SciPy calculations '
        'reproduce the R group-test p-values and BH adjustment. All coefficients, model '
        'diagnostics, standardization records, R fits and scripts are saved. The global '
        'four-hypothesis adjustment is in the [combined report](../README.md).\n'
    )
    (OUT/'README.md').write_text(report)


if __name__=='__main__':main()
