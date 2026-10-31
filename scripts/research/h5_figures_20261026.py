#!/usr/bin/env python3
"""H5 donor measurements, conventional disease box plots, and INS decomposition."""
from pathlib import Path
import os,json,itertools
os.environ.setdefault('MPLCONFIGDIR','/tmp/h5-20261026-mpl')
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('output/h5_h10_original_plan_20261026/h5')
GROUPS=['ND','Aab+','T1D'];COLORS=['#4477AA','#CC6633','#228833']

def comparisons(d,value):
    samples=[d.loc[d.donor_status.eq(g),value].dropna().to_numpy() for g in GROUPS]
    kw=stats.kruskal(*samples)
    globalrow=dict(metric=value,test='Kruskal-Wallis',statistic=kw.statistic,p_value=kw.pvalue)
    rows=[]
    for i,j in itertools.combinations(range(3),2):
        x,y=samples[i],samples[j];tied=len(np.unique(np.r_[x,y]))<len(x)+len(y)
        test=stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic' if tied else 'exact',use_continuity=True)
        rows.append(dict(metric=value,group_1=GROUPS[i],group_2=GROUPS[j],n1=len(x),n2=len(y),
                 p_value=test.pvalue,method='Tie-corrected Wilcoxon rank-sum' if tied else 'Exact Wilcoxon rank-sum'))
    q=multipletests([r['p_value'] for r in rows],method='fdr_bh')[1]
    for r,p in zip(rows,q):r['p_bh']=p
    return globalrow,rows

def boxplot(ax,d,value,ylabel,title,pairs):
    samples=[d.loc[d.donor_status.eq(g),value].dropna().to_numpy() for g in GROUPS]
    boxes=ax.boxplot(samples,positions=range(3),patch_artist=True,widths=.48,showfliers=False,
       medianprops=dict(color='#333333',lw=1.8))
    for b,c in zip(boxes['boxes'],COLORS):b.set(facecolor=c,alpha=.18,edgecolor=c)
    for i,(s,c) in enumerate(zip(samples,COLORS)):
        jitter=np.linspace(-.11,.11,len(s))
        ax.scatter(i+jitter,s,s=42,c=c,edgecolors='white',linewidths=.7,zorder=3)
    allv=np.concatenate(samples);lo=float(allv.min());hi=float(allv.max());span=max(hi-lo,abs(hi)*.15,.01)
    level=0
    for r in pairs:
        if r['p_bh']>=.05:continue
        a,b=GROUPS.index(r['group_1']),GROUPS.index(r['group_2']);y=hi+span*(.12+.16*level)
        ax.plot([a,a,b,b],[y-span*.035,y,y,y-span*.035],color='#222222',lw=1)
        p=r['p_bh'];star='***' if p<.001 else '**' if p<.01 else '*'
        ax.text((a+b)/2,y+span*.01,star,ha='center',va='bottom',fontsize=13);level+=1
    ax.set_ylim(lo-span*.12,hi+span*(.17+.16*level));ax.set_xticks(range(3),GROUPS)
    ax.set_ylabel(ylabel);ax.set_title(title,fontsize=13);ax.spines[['top','right']].set_visible(False)

def main():
    d=pd.read_parquet(OUT/'object_ins_features.parquet');rows=[]
    for (case,group),s in d.groupby(['case_id','donor_status']):
        area=s.core_area_um2.sum()/1e6
        row=dict(case_id=case,donor_status=group,endocrine_objects=len(s),
             beta_positive_objects=int(s.n_beta_core.gt(0).sum()),core_area_mm2=area,n_beta=int(s.n_beta_core.sum()),
             mean_beta_INS=s.ins_sum_beta_core.sum()/s.n_beta_ins_measured.sum(),
             beta_density=s.n_beta_core.sum()/area,
             ins_beta_density=s.ins_sum_beta_core.sum()/area,
             ins_nonbeta_density=s.ins_sum_nonbeta_core.sum()/area,
             ins_hormone_negative_density=s.ins_sum_CHGA_hormone_negative_core.sum()/area,
             ins_other_nonbeta_density=(s.ins_sum_nonbeta_core.sum()-s.ins_sum_CHGA_hormone_negative_core.sum())/area,
             ins_total_density=s.ins_sum_all_core.sum()/area,
             hormone_negative_cell_density=s.n_CHGA_hormone_negative.sum()/area,
             ins_excluded_source_density=s.ins_source_excluded_sum.sum()/area,
             ins_unresolved_source_density=s.ins_source_unresolved_margin_sum.sum()/area)
        for column in ['n_all_ins_measured','n_all_ins_missing','n_beta_ins_measured','n_beta_ins_missing',
                       'n_source_excluded_core','n_source_excluded_INS_measured','n_source_excluded_INS_missing']:
            row[column]=int(s[column].sum())
        for marker in ['NaKATPase','BActin','ECAD','BCatenin','CD44','HLADR','CHGA']:
            valid=s[f'beta_mean_{marker}'].notna()
            row[f'mean_beta_{marker}']=np.average(s.loc[valid,f'beta_mean_{marker}'],weights=s.loc[valid,'n_beta_core']) if valid.any() else np.nan
        rows.append(row)
    donor=pd.DataFrame(rows)
    donor.to_csv(OUT/'plot_data_donor_ins_measurements.csv',index=False)
    metrics=['beta_density','mean_beta_INS','ins_nonbeta_density','hormone_negative_cell_density',
             'ins_beta_density','ins_excluded_source_density']
    globaltests=[];allpairs=[]
    for metric in metrics:
        g,p=comparisons(donor,metric);globaltests.append(g);allpairs+=p
    pd.DataFrame(globaltests).to_csv(OUT/'donor_measurement_global_tests.csv',index=False)
    pd.DataFrame(allpairs).to_csv(OUT/'donor_measurement_pairwise_tests.csv',index=False)
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axs=plt.subplots(1,3,figsize=(13,4.5))
    labels=[('beta_density','β cells/mm²','β-cell abundance'),
            ('mean_beta_INS','Mean INS signal (scaled)','Retained β-cell insulin'),
            ('hormone_negative_cell_density','Cells/mm²','CHGA-positive hormone-negative cells')]
    for ax,(v,y,t) in zip(axs,labels):boxplot(ax,donor,v,y,t,[p for p in allpairs if p['metric']==v])
    fig.tight_layout()
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'insulin_quantity_and_state.{ext}',dpi=180)
    plt.close(fig)
    # All measured core signal is retained in the plotted additive decomposition.
    components=[('ins_beta_density','Typed β cells','#4477AA'),
       ('ins_hormone_negative_density','CHGA+ hormone-negative','#AA66AA'),
       ('ins_other_nonbeta_density','Other analysed cells','#AAAAAA'),
       ('ins_excluded_source_density','Source residual cells','#CCAA55')]
    ordered=pd.concat([donor[donor.donor_status.eq(g)].sort_values('case_id') for g in GROUPS])
    long=[];fig,ax=plt.subplots(figsize=(12,4.8));bottom=np.zeros(len(ordered))
    for v,label,color in components:
        ax.bar(np.arange(len(ordered)),ordered[v],bottom=bottom,color=color,label=label,width=.8)
        for (_,r),val in zip(ordered.iterrows(),ordered[v]):long.append(dict(case_id=r.case_id,donor_status=r.donor_status,component=label,signal_density=val))
        bottom+=ordered[v].to_numpy()
    ax.set_xticks(np.arange(len(ordered)),ordered.case_id.astype(str),rotation=65)
    ax.set_ylabel('INS signal / core area\n(scaled cell means/mm²)');ax.set_xlabel('Donor')
    ax.set_title('Sources of core insulin signal');ax.spines[['top','right']].set_visible(False)
    ax.legend(frameon=False,ncol=2,fontsize=9)
    for x in [6.5,13.5]:ax.axvline(x,color='#777777',lw=.7,linestyle=':')
    fig.tight_layout()
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'insulin_signal_decomposition.{ext}',dpi=180)
    plt.close(fig);pd.DataFrame(long).to_csv(OUT/'plot_data_insulin_decomposition.csv',index=False)
    slopes=pd.read_csv(OUT/'plot_data_donor_coefficients.csv');pairs=pd.read_csv(OUT/'donor_disease_pairwise_tests.csv')
    fig,ax=plt.subplots(figsize=(6.5,5.1))
    boxplot(ax,slopes,'standardized_beta','Standardized regression coefficient (β)',
            'Immune localization and β-cell insulin',pairs.to_dict('records'))
    ax.axhline(0,color='#888888',lw=1,ls='--');ax.set_xlabel('Disease state');fig.tight_layout()
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'immune_localization_beta_insulin.{ext}',dpi=190)
    plt.close(fig)
    from h5_factual_narratives_20261026 import refresh_narratives
    refresh_narratives()
    print(pd.DataFrame(globaltests).to_string(index=False));print(pd.DataFrame(allpairs).to_string(index=False))

if __name__=='__main__':main()
