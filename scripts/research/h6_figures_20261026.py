#!/usr/bin/env python3
"""All H6 null, conditional and sensitivity figures, including null results."""
from pathlib import Path
import os,itertools,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/h6-20261026-mpl')
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('output/h5_h10_original_plan_20261026/h6')
GROUPS=['ND','Aab+','T1D'];COLORS=['#4477AA','#CC6633','#228833']

def save(fig,name):
    fig.tight_layout()
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=185,bbox_inches='tight')
    plt.close(fig)

def box(ax,d,col,title,ylabel):
    samples=[d.loc[d.donor_status.eq(g),col].dropna().to_numpy() for g in GROUPS]
    bp=ax.boxplot(samples,positions=range(3),widths=.48,patch_artist=True,showfliers=False,
                  medianprops={'color':'#333333','lw':1.8})
    for i,(s,c,b) in enumerate(zip(samples,COLORS,bp['boxes'])):
        b.set(facecolor=c,edgecolor=c,alpha=.18)
        if len(s)<2:
            b.set_visible(False);bp['medians'][i].set_visible(False)
            for key in ['caps','whiskers']:
                for artist in bp[key][2*i:2*i+2]:artist.set_visible(False)
        ax.scatter(i+np.linspace(-.11,.11,len(s)),s,s=42,c=c,edgecolors='white',lw=.6,zorder=3)
    ax.set_xticks(range(3),GROUPS);ax.set_title(title);ax.set_ylabel(ylabel);ax.spines[['top','right']].set_visible(False)

def signed(v):
    v=np.asarray(v,float);v=v[np.isfinite(v)];nonzero=v[v!=0]
    if len(nonzero)<2:return np.nan
    method='exact' if len(np.unique(abs(nonzero)))==len(nonzero) else 'approx'
    return stats.wilcoxon(v,method=method,correction=True).pvalue

def paired_plot(d,left,right,name,title,ylab):
    fig,axes=plt.subplots(1,3,figsize=(11.2,4.7),sharey=True)
    for ax,g,c in zip(axes,GROUPS,COLORS):
        s=d[d.donor_status.eq(g)].dropna(subset=[left,right])
        for _,r in s.iterrows():
            ax.plot([0,1],[r[left],r[right]],color=c,alpha=.5,lw=1)
            ax.scatter([0,1],[r[left],r[right]],color=c,s=35,edgecolors='white',lw=.5,zorder=3)
        ax.set_xticks([0,1],['Observed','Permuted'] if left=='observed_percent'else['Nonpositive','Positive'])
        if left!='observed_percent':ax.set_xlabel('Excess adjacency')
        ax.set_title(g);ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel(ylab);fig.suptitle(title,y=1.02)
    save(fig,name)

def main():
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','pdf.fonttype':42})
    d=pd.read_parquet(OUT/'h6_object_features.parquet');rows=[]
    coverage=pd.read_csv(OUT/'plot_data_observed_adjacency.csv')
    for _,r in coverage.iterrows():
        s=d[d.case_id.eq(r.case_id)]
        rows.append(dict(case_id=r.case_id,donor_status=r.donor_status,objects=len(s),
            observed_percent=100*s.combined_observed.mean(),null_percent=100*s.combined_null_mean.mean(),
            excess_fraction=s.combined_excess.mean(),excess_percentage_points=100*s.combined_excess.mean(),
            alpha_excess=s.alpha_excess.mean(),delta_excess=s.delta_excess.mean()))
    donor=pd.DataFrame(rows);donor.to_csv(OUT/'plot_data_permutation_donor_summary.csv',index=False)
    tests=[]
    for g in GROUPS:
        s=donor[donor.donor_status.eq(g)].dropna(subset=['excess_fraction'])
        tests.append(dict(donor_status=g,n=len(s),test='Paired Wilcoxon signed-rank',p_value=signed(s.excess_fraction)))
    q=multipletests([r['p_value'] for r in tests],method='fdr_bh')[1]
    for r,v in zip(tests,q):r['p_bh']=v
    pd.DataFrame(tests).to_csv(OUT/'observed_vs_permuted_donor_tests.csv',index=False)
    paired_plot(donor,'observed_percent','null_percent','observed_vs_permuted_adjacency',
                'Observed and permuted endocrine adjacency','β–α/δ edges among endocrine edges (%)')
    fig,ax=plt.subplots(figsize=(6.4,5.1))
    box(ax,donor,'excess_percentage_points','Adjacency relative to random cell placement',
        'Observed minus permuted edge fraction\n(percentage points)')
    samples=[donor.loc[donor.donor_status.eq(g),'excess_percentage_points'].dropna().to_numpy() for g in GROUPS]
    kw=stats.kruskal(*samples)
    pd.DataFrame([dict(test='Kruskal-Wallis',statistic=kw.statistic,p_value=kw.pvalue,donors=19)]).to_csv(OUT/'excess_adjacency_disease_global_test.csv',index=False)
    pp=[]
    for i,j in itertools.combinations(range(3),2):
        a,b=samples[i],samples[j];tied=len(np.unique(np.r_[a,b]))<len(a)+len(b)
        p=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic'if tied else'exact',use_continuity=True).pvalue
        pp.append(dict(group_1=GROUPS[i],group_2=GROUPS[j],n1=len(a),n2=len(b),p_value=p))
    q=multipletests([r['p_value']for r in pp],method='fdr_bh')[1]
    for r,v in zip(pp,q):r['p_bh']=v
    pd.DataFrame(pp).to_csv(OUT/'excess_adjacency_disease_pairwise_tests.csv',index=False)
    lo=min(np.concatenate(samples));hi=max(np.concatenate(samples));span=hi-lo;level=0
    for r in pp:
        if r['p_bh']>=.05:continue
        a,b=GROUPS.index(r['group_1']),GROUPS.index(r['group_2']);y=hi+span*(.12+.15*level)
        ax.plot([a,a,b,b],[y-.03*span,y,y,y-.03*span],c='#222222',lw=1)
        ax.text((a+b)/2,y+.01*span,'**'if r['p_bh']<.01 else'*',ha='center',va='bottom');level+=1
    ax.set_ylim(lo-.12*span,hi+span*(.18+.15*level))
    ax.axhline(0,ls='--',c='#888888',lw=1);save(fig,'excess_adjacency_by_disease')
    # Main presentation retains independently estimated donor slopes and shows
    # the practical random-intercept model separately, without equating them.
    reg=pd.read_csv(OUT/'donor_regression_eligibility.csv')
    primary=reg[reg.model.eq('primary')&reg.estimable].copy()
    primary.to_csv(OUT/'plot_data_primary_donor_coefficients.csv',index=False)
    slopes=pd.read_csv(OUT/'disease_slopes.csv')
    key='primary_disease_interaction_intercept_only'
    selected=slopes[slopes.model.eq(key)].copy()
    selected.to_csv(OUT/'plot_data_primary_mixed_model_slopes.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(11.5,5))
    box(axes[0],primary,'standardized_beta','Independent donor estimates','Standardized regression coefficient (β)')
    for i,(g,c) in enumerate(zip(GROUPS,COLORS)):
        r=selected[selected.donor_status.eq(g)].iloc[0]
        axes[1].errorbar(i,r['contact_z.trend'],yerr=[[r['contact_z.trend']-r['lower.CL']],
                         [r['upper.CL']-r['contact_z.trend']]],fmt='o',color=c,capsize=5,ms=7)
    axes[1].set_xticks(range(3),GROUPS);axes[1].set_title('Mixed-model estimates')
    axes[1].set_ylabel('Standardized regression coefficient (β)');axes[1].spines[['top','right']].set_visible(False)
    for ax in axes:ax.axhline(0,ls='--',c='#888888',lw=1)
    fig.suptitle('Endocrine adjacency and β-cell insulin',y=1.02);save(fig,'adjacency_insulin_association')
    coeff=pd.read_csv(OUT/'model_coefficients.csv');diags=json.loads((OUT/'model_diagnostics.json').read_text())
    labels={'primary':'Combined adjacency','alpha':'Alpha adjacency','delta':'Delta adjacency',
        'beta_density':'β-cell density outcome','contact_by_immune':'Adjacency × immune burden',
        'raw_INS':'Raw INS','corrected_INS':'RedSEA-corrected INS','strong_hormone_margin':'Stronger hormone margin',
        'nerve_vessel_adjustment':'Nerve/vessel adjustment','boundary_contacts':'Cell-boundary contacts'}
    labels['nerve_vessel_proximity_adjustment']='Nerve/vessel proximity adjustment'
    sr=[]
    for name,label in labels.items():
        available=[k for k in diags if k.startswith(name+'_common_')]
        key=available[-1]
        term='contact_z:log_nearby_immune_z' if name=='contact_by_immune'else'contact_z'
        r=coeff[coeff.model.eq(key)&coeff.term.eq(term)].iloc[0].to_dict()
        sr.append(dict(analysis=name,label=label,selected_model=key,standardized_beta=r['Estimate'],lower=r['lower'],upper=r['upper'],
                       p_value=r['Pr(>|t|)'],df=r['df'],objects=diags[key]['objects'],donors=diags[key]['donors'],
                       singular=diags[key]['singular']))
    summary=pd.DataFrame(sr);summary.to_csv(OUT/'plot_data_sensitivity_estimates.csv',index=False)
    fig,ax=plt.subplots(figsize=(8.3,6.5))
    y=np.arange(len(summary))[::-1]
    ax.errorbar(summary.standardized_beta,y,xerr=[summary.standardized_beta-summary.lower,summary.upper-summary.standardized_beta],fmt='o',color='#385B75',capsize=4)
    ax.set_yticks(y,summary.label);ax.axvline(0,c='#888888',ls='--',lw=1)
    ax.set_xlabel('Standardized regression coefficient (β)');ax.set_title('Conditional associations across analysis variants')
    ax.spines[['top','right']].set_visible(False);save(fig,'adjacency_analysis_variants')
    # Prespecified coarse composition/size strata, formed without INS outcomes.
    bins={}
    for col in ['beta_fraction','delta_fraction','log_core_area']:
        d[col+'_stratum'],edges=pd.qcut(d[col],3,labels=False,duplicates='drop',retbins=True)
        bins[col]=edges.tolist()
    (OUT/'matched_strata_definitions.json').write_text(json.dumps(bins,indent=2)+'\n')
    stratumcols=['beta_fraction_stratum','delta_fraction_stratum','log_core_area_stratum']
    d['positive_excess']=d.combined_excess.gt(0)
    matched=[]
    for keys,s in d.groupby(['case_id','donor_status']+stratumcols):
        neg=s.loc[~s.positive_excess];pos=s.loc[s.positive_excess]
        if not len(neg) or not len(pos):continue
        matched.append(dict(case_id=keys[0],donor_status=keys[1],beta_stratum=keys[2],delta_stratum=keys[3],area_stratum=keys[4],
             nonpositive_objects=len(neg),positive_objects=len(pos),
             nonpositive_INS=neg.ins_mean_beta_core.mean(),positive_INS=pos.ins_mean_beta_core.mean()))
    matched=pd.DataFrame(matched);matched.to_csv(OUT/'matched_composition_size_strata.csv',index=False)
    d[['islet_id','case_id','donor_status','combined_excess','ins_mean_beta_core']+stratumcols+['positive_excess']].to_csv(OUT/'matched_strata_object_assignments.csv',index=False)
    matched_donor=matched.groupby(['case_id','donor_status'],as_index=False).agg(
        matched_strata=('beta_stratum','size'),nonpositive_objects=('nonpositive_objects','sum'),positive_objects=('positive_objects','sum'),
        nonpositive_INS=('nonpositive_INS','mean'),positive_INS=('positive_INS','mean'))
    matched_donor['INS_difference']=matched_donor.positive_INS-matched_donor.nonpositive_INS
    matched_donor.to_csv(OUT/'plot_data_matched_composition_size.csv',index=False)
    paired_plot(matched_donor,'nonpositive_INS','positive_INS','matched_composition_size_insulin',
                'Insulin within comparable composition and size strata','Mean β-cell INS signal (scaled)')
    print(summary.to_string(index=False));print(pd.DataFrame(tests).to_string(index=False));print(matched_donor.to_string(index=False))

if __name__=='__main__':main()
