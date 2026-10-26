#!/usr/bin/env python3
"""Standard scatter plots, paired donor summaries, and H4 result report."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h4-20261025-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon,spearmanr

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h1_h4_literature_methods_20261025/h4'
GROUPS=['ND','Aab+','T1D'];COLORS={'ND':'#4477AA','Aab+':'#CC6633','T1D':'#228833'}

def main():
 d=pd.read_parquet(OUT/'object_source_features.parquet')
 d=d.loc[d.variant.eq('continuous')&d.source_pool_nonzero].copy()
 r=pd.read_csv(OUT/'mixed_model_results.csv')
 primary=r.loc[r.variant.eq('continuous')&r.cohort.eq('all_complete')&r.outcome.eq('raw_density')]
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                      'svg.fonttype':'none','pdf.fonttype':42})
 fig,axs=plt.subplots(1,3,figsize=(14,4.8))
 for ax,g in zip(axs,GROUPS):
  sub=d.loc[d.donor_status.eq(g)];rr=primary.loc[primary.group.eq(g)].iloc[0]
  ax.scatter(sub.myeloid_fraction*100,sub.beta_density_mm2,s=13,alpha=.35,
             color=COLORS[g],linewidths=0,rasterized=True)
  ax.set_title(f'{g} | {len(sub):,} sections, {sub.case_id.nunique()} donors',fontweight='bold')
  ax.set_xlabel('Myeloid share of measured LGALS3 (%)')
  ax.set_xlim(-3,103);ax.set_ylim(-700,31500)
  ax.text(.04,.97,f'Adjusted change per +10 percentage points\n'
          f'{rr.estimate:+.1f} Beta cells/mm²\n'
          f'95% CI {rr.ci_low:+.1f} to {rr.ci_high:+.1f}; p={rr.p:.3f}',
          transform=ax.transAxes,va='top',fontsize=10,
          bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
 axs[0].set_ylabel('Beta cells/mm² of endocrine core')
 fig.suptitle('H4: myeloid LGALS3 allocation does not predict lower Beta density',fontsize=15,fontweight='bold')
 fig.text(.5,.01,'Each point is one endocrine section. Estimates use donor-random-intercept linear mixed models with source, immune, size, and class adjustments.',
          ha='center',fontsize=9)
 fig.tight_layout(rect=[0,.06,1,.91])
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'h4_source_scatter.{ext}',dpi=180)
 plt.close(fig)
 d['beta_present']=d.n_beta_core.gt(0)
 donors=d.groupby(['case_id','donor_status','beta_present'],observed=True).agg(
     sections=('islet_id','size'),mean_myeloid_fraction=('myeloid_fraction','mean'),
     mean_beta_density=('beta_density_mm2','mean')).reset_index()
 donors.to_csv(OUT/'donor_beta_presence_source_summary.csv',index=False)
 tests=[]
 for g in GROUPS:
  q=donors.loc[donors.donor_status.eq(g)].pivot(index='case_id',columns='beta_present',values='mean_myeloid_fraction').dropna()
  w=wilcoxon(q[True],q[False],method='exact')
  tests.append(dict(group=g,donors=len(q),mean_difference_pp=float((q[True]-q[False]).mean()*100),
                    wilcoxon_W=float(w.statistic),p=float(w.pvalue)))
 pd.DataFrame(tests).to_csv(OUT/'donor_paired_beta_presence_source_tests.csv',index=False)
 fig,axs=plt.subplots(1,3,figsize=(11,4.6),sharey=True)
 for ax,g,t in zip(axs,GROUPS,tests):
  q=donors.loc[donors.donor_status.eq(g)].pivot(index='case_id',columns='beta_present',values='mean_myeloid_fraction').dropna()
  for case,row in q.iterrows():
   ax.plot([0,1],[row[False]*100,row[True]*100],'-o',color=COLORS[g],alpha=.7,lw=1.2,ms=5)
  ax.set_xticks([0,1],['Beta absent','Beta present']);ax.set_xlim(-.2,1.2);ax.set_ylim(-2,100)
  ax.set_title(f"{g} | {len(q)} donors\nWilcoxon p={t['p']:.3f}")
 axs[0].set_ylabel('Donor mean myeloid LGALS3 share (%)')
 fig.suptitle('LGALS3 source allocation by Beta-cell presence',fontweight='bold',fontsize=15)
 fig.text(.5,.015,'Each line connects the same donor; means are over eligible sections. This paired comparison is unadjusted.',ha='center',fontsize=10)
 fig.tight_layout(rect=[0,.06,1,.91])
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'h4_donor_paired_source.{ext}',dpi=180)
 plt.close(fig)
 donorsall=d.groupby(['case_id','donor_status'],observed=True).agg(
     sections=('islet_id','size'),mean_myeloid_fraction=('myeloid_fraction','mean'),
     mean_beta_density=('beta_density_mm2','mean'),beta_count=('n_beta_core','sum'),
     core_area_um2=('core_area_um2','sum')).reset_index()
 donorsall.to_csv(OUT/'donor_source_summary.csv',index=False)
 print(primary.to_string(index=False));print(pd.DataFrame(tests).to_string(index=False))

if __name__=='__main__':main()
