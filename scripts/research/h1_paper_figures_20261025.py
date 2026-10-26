#!/usr/bin/env python3
"""Paired donor plots following Barlow Figure 3G; no resampling inference."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/islet-paper-figures')
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('output/h1_h4_literature_methods_20261025/h1')
COLORS={'ND':'#4477AA','Aab+':'#CC6633','T1D':'#228833'}


def figure_legend(metric, tests):
    if metric == 'beta_presence':
        opening = ('β-cell presence according to immune-cell location. Each point shows the percentage '
                   'of eligible endocrine objects (islets or clusters) containing at least one typed β cell '
                   'for one donor. The denominator is all eligible objects in that donor–location category.')
    else:
        opening = ('Core β-cell density according to immune-cell location. Each point shows one donor’s '
                   'arithmetic mean of individual endocrine-object β-cell densities, calculated as the '
                   'typed β-cell count divided by core area (cells/mm²). Objects without β cells contribute zero; '
                   'this is not total β-cell count divided by total tissue area.')
    pvalues = '; '.join(f"{r['group']}, p={r['p']:.6g}" for r in tests if r['metric'] == metric)
    return (opening + ' Lines connect the two measurements from the same donor. “Peri only” indicates '
            'at least one immune cell in the surrounding 20-µm band and none in the endocrine core. '
            '“Core present” indicates at least one core immune cell, whether or not peri immune cells '
            'are also present. The peri band excludes other endocrine cores and is clipped to tissue. '
            'Only immune-exposed objects passing the corrected counting and acquisition-edge exclusions '
            'are included: ND, 792 objects from 7 donors; Aab+, 1,057 objects from 7 donors; T1D, '
            '1,088 objects from 6 donors. ND, non-diabetic; Aab+, autoantibody-positive; T1D, type 1 diabetes. '
            'Two-sided paired Wilcoxon signed-rank tests compare immune-location categories within each '
            f'disease state: {pvalues}. P values are unadjusted. *p<0.05; no asterisk denotes p≥0.05. '
            'These paired comparisons are supplementary to the direct comparisons across disease states.')


def main():
    d=pd.read_csv(OUT/'analysis_objects.csv',dtype={'case_id':str})
    d['beta_present']=d.n_beta_core.gt(0)
    d['core_immune']=d.n_immune_core.gt(0)
    donor=d.groupby(['donor_status','case_id','core_immune']).agg(
        objects=('islet_id','size'),beta_density=('beta_density_mm2','mean'),
        beta_presence=('beta_present','mean')).reset_index()
    donor.to_csv(OUT/'paired_donor_measurements.csv',index=False)
    tests=[]
    legends=[]
    for metric in ['beta_density','beta_presence']:
        fig,axs=plt.subplots(1,3,figsize=(10.8,4.6),sharey=True)
        factor=100 if metric=='beta_presence' else 1
        upper=100 if metric=='beta_presence' else float(donor[metric].max())
        bracket_y=upper*1.03
        for ax,(group,color) in zip(axs,COLORS.items()):
            x=donor.loc[donor.donor_status.eq(group)].pivot(index='case_id',columns='core_immune',values=metric).dropna()
            for i,(case,r) in enumerate(x.iterrows()):
                ax.plot([0,1],r.to_numpy()*factor,'o-',color=color,alpha=.7,lw=1.2,ms=5)
            test=stats.wilcoxon(x[True],x[False],alternative='two-sided',method='auto')
            tests.append(dict(group=group,metric=metric,donors=len(x),
                mean_core=x[True].mean(),mean_peri_only=x[False].mean(),
                mean_difference=(x[True]-x[False]).mean(),p=test.pvalue))
            ax.set_title(group,fontsize=13)
            ax.set_xticks([0,1],['Peri only','Core present'])
            ax.set_xlabel('Immune-cell location')
            ax.set_xlim(-.2,1.2);ax.set_ylim(0,upper*1.16)
            if metric=='beta_presence':ax.set_yticks(np.arange(0,101,20))
            ax.spines[['top','right']].set_visible(False)
            if test.pvalue < .05:
                ax.plot([0,0,1,1],[bracket_y-upper*.02,bracket_y,bracket_y,bracket_y-upper*.02],
                        color='black',lw=1)
                ax.text(.5,bracket_y+upper*.005,'*',ha='center',va='bottom',fontsize=18)
        axs[0].set_ylabel('Endocrine objects with ≥1 β cell (%)' if metric=='beta_presence' else
                          'Donor mean core β-cell density\n(cells/mm²)',fontsize=11)
        fig.suptitle('β-cell presence' if metric=='beta_presence' else 'Core β-cell density',fontsize=15)
        fig.tight_layout(rect=(0,0,1,.93))
        fig.savefig(OUT/f'h1_{metric}_paired.png',dpi=180)
        fig.savefig(OUT/f'h1_{metric}_paired.pdf')
        fig.savefig(OUT/f'h1_{metric}_paired.svg')
        plt.close(fig)
        legend=figure_legend(metric,tests)
        (OUT/f'h1_{metric}_paired_legend.txt').write_text(legend+'\n')
        legends.append(f'![{metric.replace("_"," ")}]({(OUT/f"h1_{metric}_paired.png").resolve()})\n\n{legend}\n')
    pd.DataFrame(tests).to_csv(OUT/'paired_donor_tests.csv',index=False)
    (OUT/'figure_legends.md').write_text('# Paired donor comparisons\n\n'+'\n'.join(legends))
    print(pd.DataFrame(tests).to_string(index=False))


if __name__=='__main__':main()
