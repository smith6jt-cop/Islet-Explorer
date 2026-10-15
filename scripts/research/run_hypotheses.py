#!/usr/bin/env python3
"""Donor-aware primary analyses and explicit unavailability; never use DPT labels.

Continuous endpoints use donor-balanced WLS with small-sample donor-clustered
inference. Proportions with integer denominators use beta-binomial likelihood
and donor bootstrap. Ten BH slots are retained when some tests are unavailable.
"""
from __future__ import annotations
import argparse
import json
import shutil
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, betaln, gammaln
from scipy.stats import norm
import statsmodels.api as sm
from common import MARKERS, SEED, sha256, write_json, require_cohort, offline_report_root

def within_between(d, col, prefix):
    d=d.copy();d[prefix+'_between']=d.groupby('case_id')[col].transform('mean')
    d[prefix+'_within']=d[col]-d[prefix+'_between'];return d

def design(d, interaction=False):
    d=within_between(d,'log_area','log_area')
    x=pd.DataFrame({'Intercept':np.ones(len(d)),'disease':(d.donor_status=='T1D').astype(float),
                    'log_area_within':d.log_area_within,'log_area_between':d.log_area_between,
                    'age':d.age,'sex':d.sex.astype(str).str.lower().str.startswith('m').astype(float)},index=d.index)
    if (d.donor_status=='Aab+').any():
        x['aab_status']=(d.donor_status=='Aab+').astype(float)
    if interaction:
        x['disease:log_area_within']=x.disease*x.log_area_within
        x['disease:log_area_between']=x.disease*x.log_area_between
    return x

def beta_binomial_fit(x,y,n,weights=None):
    X=np.asarray(x,float);y=np.asarray(y,float);n=np.asarray(n,float)
    if np.any(y<0)|np.any(y>n)|np.any(n<=0): raise ValueError('Invalid binomial counts')
    weights=np.ones(len(y)) if weights is None else weights
    # Standardize nonconstant predictors for numerical conditioning; map back.
    mu=X.mean(axis=0);sd=X.std(axis=0);mu[sd==0]=0;sd[sd==0]=1
    Z=(X-mu)/sd
    def objective(theta):
        p=np.clip(expit(Z@theta[:-1]),1e-9,1-1e-9);precision=np.exp(theta[-1])
        a=p*precision;b=(1-p)*precision
        ll=betaln(y+a,n-y+b)-betaln(a,b)+gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1)
        return -np.sum(weights*ll)
    result=minimize(objective,np.r_[np.zeros(X.shape[1]),np.log(10)],method='L-BFGS-B',bounds=[(None,None)]*X.shape[1]+[(-6,16)])
    if not result.success or not np.isfinite(result.fun): raise ValueError('Beta-binomial fit did not converge')
    beta=result.x[:-1]/sd;beta[0]-=np.sum(mu*beta)
    return beta

def donor_bootstrap(d, fitter, repetitions=200):
    rng=np.random.default_rng(SEED);donors=np.array(sorted(d.case_id.unique()));values=[]
    for _ in range(repetitions):
        parts=[]
        for slot,cid in enumerate(rng.choice(donors,len(donors),replace=True)):
            sub=d.loc[d.case_id==cid].copy();sub['case_id']=f'{cid}@{slot}';parts.append(sub)
        try:
            value=fitter(pd.concat(parts,ignore_index=True))
            if np.isfinite(value):values.append(float(value))
        except (ValueError,np.linalg.LinAlgError):pass
    return values

def continuous(d,columns,term,bootstrap=200):
    def fit(sub):
        X=design(sub)
        for name in columns:
            centered=within_between(sub,name,name)
            X[name+'_within']=centered[name+'_within'].to_numpy()
            X[name+'_between']=centered[name+'_between'].to_numpy()
        if np.linalg.matrix_rank(X)<X.shape[1]: raise ValueError('Rank-deficient design')
        w=1/sub.groupby('case_id').case_id.transform('size')
        return sm.WLS(sub.response,X,weights=w).fit(cov_type='cluster',cov_kwds={'groups':sub.case_id,'use_correction':True},use_t=True)
    result=fit(d); draws=donor_bootstrap(d,lambda x:fit(x).params[term],bootstrap)
    if len(draws)<.8*bootstrap: raise ValueError('Fewer than 80% donor bootstrap fits converged')
    ci=np.quantile(draws,[.025,.975]);return dict(effect=float(result.params[term]),ci_low=float(ci[0]),ci_high=float(ci[1]),p_value=float(result.pvalues[term]),bootstrap_success=len(draws))

def proportion(d,bootstrap=200,exposure=None,extra=(),size_modifier=None):
    term='disease:log_area_within' if exposure is None else exposure+'_within'
    def fit(sub):
        X=design(sub,interaction=exposure is None)
        if exposure is not None:
            sub=within_between(sub,exposure,exposure)
            X[exposure+'_within']=sub[exposure+'_within'];X[exposure+'_between']=sub[exposure+'_between']
        for col in extra:
            sub=within_between(sub,col,col)
            X[col+'_within']=sub[col+'_within'];X[col+'_between']=sub[col+'_between']
        if size_modifier is not None:
            X['size:structure']=X.log_area_within*sub[size_modifier]
        if np.linalg.matrix_rank(X)<X.shape[1]: raise ValueError('Rank-deficient design')
        weights=1/sub.groupby('case_id').case_id.transform('size').to_numpy();weights*=len(weights)/weights.sum()
        return beta_binomial_fit(X,sub.success,sub.trials,weights)[list(X).index('size:structure' if size_modifier is not None else term)]
    effect=float(fit(d));draws=donor_bootstrap(d,fit,bootstrap)
    if len(draws)<.8*bootstrap: raise ValueError('Fewer than 80% donor bootstrap fits converged')
    se=np.std(draws,ddof=1);ci=np.quantile(draws,[.025,.975])
    return dict(effect=effect,ci_low=float(ci[0]),ci_high=float(ci[1]),p_value=float(2*norm.sf(abs(effect/se))) if se>0 else 1.,bootstrap_success=len(draws))

def bh_ten(pvalues):
    p=np.array(pvalues,float);ok=np.isfinite(p);q=np.full(len(p),np.nan)
    order=np.flatnonzero(ok)[np.argsort(p[ok])]
    adjusted=np.minimum.accumulate((p[order]*10/np.arange(1,len(order)+1))[::-1])[::-1]
    q[order]=np.minimum(1,adjusted);return q

def run_primary(root,registry,bootstrap,report_root):
    regions=pd.read_parquet(root/'analysis/regions.parquet')
    donors=pd.read_csv(root/'donor_metadata.csv',dtype={'case_id':str})
    require_cohort(regions,donors.case_id)
    comp=[]
    for cid in donors.case_id:
        d=pd.read_csv(root/'regional'/f'islet_composition_{int(cid):04d}.tsv',sep='\t')
        d['islet_id']=cid+'_'+d['name'];comp.append(d)
    comp=pd.concat(comp,ignore_index=True)
    comp['success']=comp.Ins_any;comp['trials']=comp.cells_total-comp.Triple_Neg
    if (comp.success>comp.trials).any():raise ValueError('Endocrine union smaller than INS-positive count')
    result=[];endpoint_tables={};sensitivity=[]
    reasons={'H4':'Validated CD45-positive cell gate and PDPN lymphatic validation unavailable.',
             'H6':'CD3e/CD8a parent-chain gate validation unavailable.',
             'H7':'No validated ColIV boundary mask or CD45 gate.',
             'H8':'No validated CK19/panCK duct mask; negative-only Ductal rule is not specific.',
             'H9':'Endocrine-cell-specific 31-marker profile and segmentation/technical dispersion validation pending.',
             'H10':'Complete validated spatial feature set and endocrine composition baseline unavailable; nested classifier must not silently omit them.'}
    def endpoint(h,band=20,flagged=False,exclude_hpap=False):
        r=regions.loc[(regions.band_um==band)&regions.geometry_valid&(regions.valid_area_um2>0)&(regions.core_tissue_coverage>.99)].copy()
        if not flagged:r=r.loc[~r.overlap&(r.annotation_review!='pending')]
        r=r.loc[r.donor_status.isin(['ND','T1D'])]
        if exclude_hpap:r=r.loc[r.donor_source!='HPAP']
        core=r.loc[r.compartment=='core'].copy();peri=r.loc[r.compartment=='peri'].copy()
        if h=='H1': d=core.merge(comp[['islet_id','success','trials']],on='islet_id',validate='one_to_one');d=d.loc[d.trials>0]
        elif h in ['H2','H3']:
            field='nerve_area_fraction' if h=='H2' else 'capillary_area_fraction'
            d=core.merge(peri[['islet_id',field]],on='islet_id',suffixes=('_core','_peri'),validate='one_to_one');d['response']=d[field+'_peri']-d[field+'_core']
        else:raise ValueError('No certified endpoint')
        d=d.merge(donors[['case_id','age','sex']],on='case_id',validate='many_to_one');d['log_area']=np.log(d.area_um2)
        d=d.dropna(subset=['age','sex','log_area'])
        return d
    for h in registry['hypotheses']:
        rec=dict(id=h['id'],title=h['title'],status='unavailable',reason='',effect=np.nan,ci_low=np.nan,ci_high=np.nan,p_value=np.nan,
                 n_donors=0,n_islets=0,unit='log odds interaction' if h['id']=='H1' else 'area fraction interaction',bootstrap_success=0)
        if h['id'] in ['H1','H2','H3']:
            d=endpoint(h['id']);endpoint_tables[h['id']]=d
            fit=lambda df,b=bootstrap: proportion(df,b) if h['id']=='H1' else continuous(df,[],'disease',b)
            try:
                rec.update(fit(d));rec.update(status='estimated_pending_scientific_review',n_donors=int(d.case_id.nunique()),n_islets=int(d.islet_id.nunique()))
                for name,band,flagged,hpap in [('band10',10,False,False),('band30',30,False,False),('flagged',20,True,False),('exclude_HPap',20,False,True)]:
                    if h['id']=='H1' and band!=20:continue
                    try:
                        s=endpoint(h['id'],band,flagged,hpap);r=fit(s,bootstrap)
                        sensitivity.append(dict(id=h['id'],analysis=name,n_donors=s.case_id.nunique(),n_islets=len(s),**r))
                    except (ValueError,np.linalg.LinAlgError) as e:sensitivity.append(dict(id=h['id'],analysis=name,status=str(e)))
                for donor in sorted(d.case_id.unique()):
                    sub=d.loc[d.case_id!=donor]
                    try:
                        if h['id']=='H1':
                            X=design(sub,True);beta=beta_binomial_fit(X,sub.success,sub.trials)
                            value=beta[list(X).index('disease:log_area_within')]
                        else:
                            X=design(sub);value=sm.WLS(sub.response,X,weights=1/sub.groupby('case_id').case_id.transform('size')).fit().params['disease']
                        sensitivity.append(dict(id=h['id'],analysis='omit_'+donor,effect=value,n_donors=sub.case_id.nunique(),n_islets=len(sub)))
                    except (ValueError,np.linalg.LinAlgError):pass
            except (ValueError,np.linalg.LinAlgError) as e:rec['reason']=str(e)
        elif h['id']=='H5':
            rec['reason']='Analytic uniform-location null is available as descriptive supporting output; required shape-preserving spatial-null validation has not been established.'
        else:rec['reason']=reasons[h['id']]
        result.append(rec);print(h['id'],rec['status'],rec['p_value'],flush=True)
    results=pd.DataFrame(result);results['q_value']=bh_ten(results.p_value);results['primary_family']='ten_hypotheses_v1'
    for h,d in endpoint_tables.items():d.to_csv(report_root/f'{h}-endpoint.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(report_root/'sensitivities.csv',index=False)
    return results,endpoint_tables

def followups(results, tables=None, bootstrap=200):
    r=results.set_index('id');estimated={i for i in r.index if np.isfinite(r.loc[i,'q_value'])};supported={i for i in estimated if r.loc[i,'q_value']<.05}
    rows=[]
    if {'H1','H2','H3'}<=estimated:
        structure=min(['H2','H3'],key=lambda i:(r.loc[i,'q_value'],0 if i=='H3' else 1))
        branch='size_by_'+('nerve' if structure=='H2' else 'capillary') if 'H1' in supported and ({'H2','H3'}&supported) else 'H1_adjust_both_structures'
        rec=dict(id='F1',branch=branch,status='exploratory_same_cohort')
        if tables is not None:
            d=tables['H1'].copy()
            try:
                rec.update(proportion(d,bootstrap,extra=['nerve_area_fraction','capillary_area_fraction'],
                                      size_modifier=('nerve_area_fraction' if structure=='H2' else 'capillary_area_fraction') if branch.startswith('size_by_') else None))
            except ValueError as e:rec.update(status='unavailable',reason=str(e))
        rows.append(rec)
    else:rows.append(dict(id='F1',branch='',status='awaiting required primary results'))
    for i,needed in [('F2',{'H4'}),('F3',{'H7','H8'}),('F4',{'H9','H10'})]:
        rows.append(dict(id=i,branch='',status='awaiting required primary results' if not needed<=estimated else 'exploratory follow-up pending'))
    return pd.DataFrame(rows)

def plots(root,tables):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for h,d in tables.items():
        d=d.copy();d['response']=d.success/d.trials if h=='H1' else d.response
        fig,ax=plt.subplots(figsize=(7,4))
        for status,g in d.groupby('donor_status'):
            ax.scatter(np.log(g.area_um2),g.response,s=4,alpha=.25,label=status)
        ax.set(xlabel='log section area (um²)',ylabel='INS-positive / endocrine-positive union' if h=='H1' else 'Peri minus core area fraction',title=h+' — descriptive islet measurements')
        ax.legend();fig.tight_layout()
        for suffix in ['png','pdf']:fig.savefig(root/f'{h}.{suffix}',dpi=150)
        plt.close(fig)
        examples=d.assign(distance=(d.response-d.response.median()).abs()).sort_values(['distance','islet_id']).groupby('donor_status').head(1)
        examples[['islet_id','case_id','islet_key','donor_status']].to_csv(root/f'{h}-representative-islets.csv',index=False)

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--report-root',required=True);p.add_argument('--registry',required=True);p.add_argument('--bootstrap',type=int,default=200);a=p.parse_args()
    root=Path(a.input_root);out=offline_report_root(root,a.report_root);registry=json.loads(Path(a.registry).read_text())
    if (out/'results.csv').exists():raise ValueError('Use a fresh primary report directory; preserve existing results')
    shutil.copyfile(a.registry,out/'registry.json')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        results,tables=run_primary(root,registry,a.bootstrap,out)
    results.to_csv(out/'results.csv',index=False);followups(results,tables,a.bootstrap).to_csv(out/'followups.csv',index=False);plots(out,tables)
    write_json(out/'run.json',dict(registry_sha256=sha256(a.registry),cohort_sha256=sha256(root/'cohort.csv'),
               runner_sha256=sha256(__file__),bootstrap_requested=a.bootstrap,completed_tests=int(results.p_value.notna().sum()),
               primary_family_size=10,acquisition='One acquisition per donor; acquisition and donor effects are inseparable.',
               scientific_acceptance=False,pseudotime_used=False))

if __name__=='__main__':main()
