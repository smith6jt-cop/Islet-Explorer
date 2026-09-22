#!/usr/bin/env python3
"""Reanalysis of canonical post-typing endpoints; no legacy export fallback.

All data products are regenerated from the selected geometry/typing revision.
One primary test per registered hypothesis, with unavailable prerequisites
retaining their BH slot. Reports are written outside the application data root.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import norm
from joblib import Parallel,delayed
from threadpoolctl import threadpool_limits
from common import MARKERS,SEED,sha256,write_json,offline_report_root
from composition import canonical_cells
from run_hypotheses import design,within_between,beta_binomial_fit,bh_ten
from prediction import permutation_test


def model_matrix(d,exposure=None,interaction=False,extra=(),product=None,disease_interaction=False):
    X=design(d,interaction=interaction)
    for name in ([exposure] if exposure else [])+list(extra):
        z=within_between(d,name,name)
        X[name+'_within']=z[name+'_within'];X[name+'_between']=z[name+'_between']
    if product:
        a,b=product
        X['product']=X[a]*X[b]
    if disease_interaction:
        X['disease:exposure']=X.disease*X[exposure+'_within']
    if not np.isfinite(X.to_numpy(float)).all() or np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError('Nonfinite or rank-deficient donor-aware design')
    return X


def negative_binomial_fit(X,y,offset,weights):
    X=np.asarray(X,float);y=np.asarray(y,float);offset=np.asarray(offset,float)
    if np.any(y<0) or np.any(y!=np.floor(y)) or not np.isfinite(offset).all():raise ValueError('Invalid count or area offset')
    center=X.mean(axis=0);scale=X.std(axis=0);center[scale==0]=0;scale[scale==0]=1;Z=(X-center)/scale
    def objective(theta):
        mu=np.exp(np.clip(Z@theta[:-1]+offset,-30,30));r=np.exp(-theta[-1])
        ll=gammaln(y+r)-gammaln(r)-gammaln(y+1)+r*np.log(r/(r+mu))+y*np.log(mu/(r+mu))
        return -float(np.sum(weights*ll))
    init=np.zeros(X.shape[1]+1);init[0]=np.log((y.sum()+.5)/np.exp(offset).sum());init[-1]=0
    result=minimize(objective,init,method='L-BFGS-B',bounds=[(None,None)]*X.shape[1]+[(-12,12)])
    if not result.success or not np.isfinite(result.fun):raise ValueError('Negative-binomial fit did not converge')
    beta=result.x[:-1]/scale;beta[0]-=np.sum(center*beta);return beta


def estimate(d,family,term,exposure=None,interaction=False,extra=(),product=None,disease_interaction=False):
    X=model_matrix(d,exposure,interaction,extra,product,disease_interaction)
    if term not in X:raise ValueError(f'Model lacks primary term {term}')
    weights=1/d.groupby('case_id').case_id.transform('size').to_numpy();weights*=len(weights)/weights.sum()
    if family=='beta-binomial':beta=beta_binomial_fit(X,d.success,d.trials,weights)
    elif family=='negative-binomial':beta=negative_binomial_fit(X,d.response,np.log(d.offset_area),weights)
    elif family=='continuous':beta=np.asarray(sm.WLS(d.response,X,weights=weights).fit().params)
    else:raise ValueError(f'Unknown endpoint family {family}')
    return float(beta[list(X).index(term)])


def bootstrap_model(d,spec,repetitions=999,workers=1,checkpoint_root=None):
    effect=estimate(d,**spec);donors=np.array(sorted(d.case_id.unique()));rng=np.random.default_rng(SEED)
    draws=[rng.choice(donors,len(donors),replace=True) for _ in range(repetitions)]
    checkpoint=Path(checkpoint_root) if checkpoint_root else None
    fingerprint=hashlib.sha256(pd.util.hash_pandas_object(d,index=False).to_numpy().tobytes()+json.dumps(spec,sort_keys=True).encode()+sha256(__file__).encode()).hexdigest()
    if checkpoint:checkpoint.mkdir(parents=True,exist_ok=True)
    def one(index,draw):
        path=checkpoint/f'{index}.json' if checkpoint else None
        if path and path.exists():
            saved=json.loads(path.read_text())
            if saved.get('fingerprint')==fingerprint:return np.nan if saved['value'] is None else saved['value']
        parts=[]
        for j,cid in enumerate(draw):
            part=d.loc[d.case_id==cid].copy();part.case_id=f'{cid}@{j}';parts.append(part)
        try:
            with threadpool_limits(limits=1):value=estimate(pd.concat(parts,ignore_index=True),**spec)
        except (ValueError,np.linalg.LinAlgError):value=np.nan
        if path:write_json(path,dict(fingerprint=fingerprint,value=float(value) if np.isfinite(value) else None))
        return value
    values=np.array(Parallel(n_jobs=workers)(delayed(one)(index,draw) for index,draw in enumerate(draws)));good=values[np.isfinite(values)]
    if len(good)<max(2,.8*repetitions):raise ValueError(f'Only {len(good)}/{repetitions} donor bootstrap fits converged')
    lo,hi=np.quantile(good,[.025,.975]);se=np.std(good,ddof=1)
    return dict(effect=effect,ci_low=float(lo),ci_high=float(hi),p_value=float(2*norm.sf(abs(effect/se))) if se>0 else 1.,
                bootstrap_requested=repetitions,bootstrap_success=len(good),bootstrap_failed=repetitions-len(good))


def endocrine_profiles(summaries):
    s=summaries.loc[summaries.compartment.eq('core')&summaries.population.eq('Endocrine union')].copy()
    if s.duplicated(['islet_id','marker']).any():raise ValueError('Duplicate endocrine marker summaries')
    profiles=s.pivot(index=['islet_id','case_id'],columns='marker',values='mean_log1p').reset_index()
    n=s.groupby('islet_id').n_cells.min();profiles['endocrine_measured_cells']=profiles.islet_id.map(n)
    profiles['profile_eligible']=(profiles.endocrine_measured_cells>=5)&np.isfinite(profiles[MARKERS]).all(axis=1)
    return profiles


def joint_counts(root,memberships,donors):
    rows=[]
    for cid in donors:
        cells=canonical_cells(pd.read_parquet(root/'typing/cells'/f'{cid}.parquet'),pd.read_parquet(root/'corrected'/f'{cid}-cells.parquet'))
        m=memberships.loc[memberships.case_id.eq(cid)].copy();m.object_id=m.object_id.astype(str)
        d=m.merge(cells,on=['case_id','object_id'],validate='many_to_one')
        for (iid,comp,band),g in d.groupby(['islet_id','compartment','band_um']):
            measured=g.loc[g.measured];valid=len(measured)>0 and all('positive__'+x in measured and measured['positive__'+x].notna().all() for x in ['CD3e','CD8a'])
            cd45_valid=len(measured)>0 and 'positive__CD45' in measured and measured.positive__CD45.notna().all()
            rows.append(dict(islet_id=iid,case_id=cid,compartment=comp,band_um=band,
                             cd8_count=float((measured.positive__CD3e&measured.positive__CD8a).sum()) if valid else np.nan,
                             cd45_count=float(measured.positive__CD45.sum()) if cd45_valid else np.nan))
    return pd.DataFrame(rows)


def prepare(root):
    islets=pd.read_parquet(root/'analysis/islets.parquet');islets.case_id=islets.case_id.astype(str)
    donors=pd.read_csv(root/'donor_metadata.csv',dtype={'case_id':str})
    if set(islets.case_id)!=set(donors.case_id):raise ValueError('Islet and donor cohorts do not match')
    base=islets.merge(donors[[c for c in ['case_id','age','sex','donor_status','donor_source'] if c=='case_id' or c not in islets]],on='case_id',validate='many_to_one')
    base['log_area']=np.log(base.area_um2)
    summary=pd.read_parquet(root/'analysis/marker_summaries.parquet');summary.case_id=summary.case_id.astype(str)
    counts=pd.read_parquet(root/'trajectory/endocrine_counts.parquet')
    base=base.merge(counts[['islet_id','endocrine_count']+[c for c in counts if c.startswith('endocrine_component__')]],on='islet_id',validate='one_to_one')
    ins=summary.loc[summary.compartment.eq('core')&summary.population.eq('Endocrine union')&summary.marker.eq('INS')]
    base=base.merge(ins[['islet_id','positive_count','positive_available']].rename(columns={'positive_count':'success','positive_available':'ins_available'}),on='islet_id',validate='one_to_one')
    base['trials']=base.endocrine_count
    if (base.success>base.trials).any():raise ValueError('INS-positive endocrine count exceeds its expanded union')
    regions=pd.read_parquet(root/'analysis/regions.parquet');regions.case_id=regions.case_id.astype(str)
    return base,regions,summary,endocrine_profiles(summary),donors


def region_wide(regions,band=20,flagged=False):
    r=regions.loc[regions.band_um.eq(band)&regions.geometry_valid&(regions.valid_area_um2>0)].copy()
    if not flagged and 'overlap' in r:r=r.loc[~r.overlap]
    if not flagged and 'primary_assessable' in r:r=r.loc[r.primary_assessable.fillna(False)]
    if 'annotation_status' in r:r=r.loc[r.annotation_status.ne('technical_region_unassessable')]
    if 'core_tissue_coverage' in r:r=r.loc[r.core_tissue_coverage>.99]
    # A missing mask is never a biological zero. Each positive-area source must
    # carry affirmative validation, independently of whether its value is zero.
    for prefix in ['nerve','capillary','lymphatic']:
        value=prefix+'_area_fraction';flag=prefix+'_validated'
        if value not in r:r[value]=np.nan
        r.loc[~r.get(flag,pd.Series(False,index=r.index)).fillna(False),value]=np.nan
    for value,flag in [('coliv_boundary_coverage','coliv_validated'),('duct_distance_um','duct_validated')]:
        if value not in r:r[value]=np.nan
        r.loc[~r.get(flag,pd.Series(False,index=r.index)).fillna(False),value]=np.nan
    fields=['valid_area_um2','nerve_area_fraction','capillary_area_fraction','lymphatic_area_fraction','coliv_boundary_coverage','duct_distance_um']
    if r.duplicated(['islet_id','compartment']).any():raise ValueError('Duplicate analysis regions')
    r=r.loc[r.islet_id.isin(r.loc[r.compartment.eq('core'),'islet_id'])]
    w=r.pivot(index='islet_id',columns='compartment',values=fields).reindex(columns=pd.MultiIndex.from_product([fields,['core','peri']]))
    w.columns=['_'.join(x) for x in w.columns]
    return w.reset_index()


def positive_counts(summary,marker,population='All cells'):
    s=summary.loc[summary.marker.eq(marker)&summary.population.eq(population)].copy()
    s.loc[~s.positive_available,'positive_count']=np.nan
    w=s.pivot(index='islet_id',columns='compartment',values='positive_count');w.columns=[marker+'_'+c for c in w]
    return w.reset_index()


def band_positive_counts(summary,marker,root,band):
    if band==20:return positive_counts(summary,marker)
    if marker!='CD45':raise ValueError('Requested nonprimary band marker counts have not been materialized')
    calls=pd.read_parquet(root/'analysis/joint_marker_counts.parquet')
    calls=calls.loc[(calls.compartment.eq('core')&calls.band_um.eq(20))|(calls.compartment.eq('peri')&calls.band_um.eq(band))]
    table=calls.pivot(index='islet_id',columns='compartment',values='cd45_count');table.columns=['CD45_'+c for c in table]
    return table.reset_index()


def weighted_coupling(coupling,band=20):
    required={'null_method','null_successes','nerve_area_um2','excess_proximity','compartment','band_um'}
    if not required.issubset(coupling) or not coupling.null_method.eq('shape_preserving_compartment_rigid').all():
        raise ValueError('Actual nerve area and shape-preserving compartment null required')
    rows=[]
    for iid,g in coupling.loc[coupling.band_um.eq(band)].groupby('islet_id'):
        if len(g)!=2 or set(g.compartment)!={'core','peri'}:continue
        known=g.nerve_area_um2.notna()&((g.nerve_area_um2==0)|((g.null_successes>=999)&g.excess_proximity.notna()))
        weight=g.nerve_area_um2.to_numpy(float)
        if not known.all() or weight.sum()<=0:continue
        value=float(np.sum(weight*g.excess_proximity.fillna(0).to_numpy(float))/weight.sum())
        rows.append(dict(islet_id=iid,excess_proximity=value,assessed_nerve_area_um2=float(weight.sum())))
    return pd.DataFrame(rows,columns=['islet_id','excess_proximity','assessed_nerve_area_um2'])


def endpoint(h,base,regions,summary,profiles,root,band=20,flagged=False):
    d=base.merge(region_wide(regions,band,flagged),on='islet_id',validate='one_to_one')
    if not flagged and 'primary_assessable' in d:d=d.loc[d.primary_assessable.fillna(False)]
    if 'annotation_status' in d:d=d.loc[d.annotation_status.ne('technical_region_unassessable')]
    if h not in ['H4','H6','H7','H8']:d=d.loc[d.donor_status.isin(['ND','T1D'])].copy()
    needed=['age','sex','log_area'];spec={};reason=None
    if h=='H1':spec=dict(family='beta-binomial',term='disease:log_area_within',interaction=True);needed+=['success','trials']
    elif h in ['H2','H3']:
        field='nerve' if h=='H2' else 'capillary';d['response']=d[field+'_area_fraction_peri']-d[field+'_area_fraction_core'];needed+=['response']
        spec=dict(family='continuous',term='disease')
    elif h=='H4':
        d=d.merge(band_positive_counts(summary,'CD45',root,band),on='islet_id',validate='one_to_one')
        d['response']=d.CD45_peri;d['offset_area']=d.valid_area_um2_peri;d['lymphatic']=d.lymphatic_area_fraction_peri
        needed+=['response','offset_area','lymphatic'];spec=dict(family='negative-binomial',term='lymphatic_within',exposure='lymphatic')
    elif h=='H5':
        coupling=weighted_coupling(pd.read_parquet(root/'analysis/coupling.parquet'),band)
        d=d.merge(coupling,on='islet_id',validate='one_to_one');d['response']=d.excess_proximity
        valid_area=d.valid_area_um2_core+d.valid_area_um2_peri
        d['nerve']=(d.nerve_area_fraction_core*d.valid_area_um2_core+d.nerve_area_fraction_peri*d.valid_area_um2_peri)/valid_area
        d['capillary']=(d.capillary_area_fraction_core*d.valid_area_um2_core+d.capillary_area_fraction_peri*d.valid_area_um2_peri)/valid_area
        needed+=['response','nerve','capillary']
        spec=dict(family='continuous',term='disease',extra=['nerve','capillary'])
    elif h=='H6':
        jc=pd.read_parquet(root/'analysis/joint_marker_counts.parquet');jc=jc.loc[jc.compartment.eq('core')&jc.band_um.eq(20)]
        d=d.merge(jc[['islet_id','cd8_count']],on='islet_id',validate='one_to_one');d['cd8']=d.cd8_count/d.valid_area_um2_core
        needed+=['success','trials','cd8'];spec=dict(family='beta-binomial',term='cd8_within',exposure='cd8')
    elif h=='H7':
        d=d.merge(band_positive_counts(summary,'CD45',root,band),on='islet_id',validate='one_to_one')
        d['success']=d.CD45_core;d['trials']=d.CD45_core+d.CD45_peri;d['coliv']=d.coliv_boundary_coverage_core
        needed+=['success','trials','coliv'];spec=dict(family='beta-binomial',term='coliv_within',exposure='coliv')
    elif h=='H8':
        d['duct_distance']=d.duct_distance_um_core;needed+=['success','trials','duct_distance']
        spec=dict(family='beta-binomial',term='duct_distance_within',exposure='duct_distance')
    elif h=='H9':
        d=d.merge(profiles.loc[profiles.profile_eligible],on=['islet_id','case_id'],validate='one_to_one');rows=[]
        for cid,g in d.groupby('case_id'):
            if len(g)<3:continue
            # Adjust section size within each donor. All markers contribute to
            # dispersion; degrees-of-freedom correct the residual variance.
            X=np.c_[np.ones(len(g)),g.log_area-g.log_area.mean()];Y=g[MARKERS].to_numpy(float)
            residual=Y-X@np.linalg.lstsq(X,Y,rcond=None)[0]
            rec=g.iloc[0][['case_id','donor_status','age','sex']].to_dict()
            variance=np.sum(residual**2,axis=0)/(len(g)-np.linalg.matrix_rank(X))
            rec.update(response=float(np.sqrt(variance.mean())),n_islets=len(g),log_area=float(g.log_area.mean()))
            rec.update({'dispersion__'+m:float(v) for m,v in zip(MARKERS,np.sqrt(variance))})
            rows.append(rec)
        return pd.DataFrame(rows),dict(family='donor-dispersion',term='disease')
    else:raise ValueError(f'Unknown direct endpoint {h}')
    d=d.dropna(subset=needed)
    if spec['family']=='beta-binomial':d=d.loc[(d.trials>0)&(d.success>=0)&(d.success<=d.trials)]
    if spec['family']=='negative-binomial':d=d.loc[d.offset_area>0]
    if d.case_id.nunique()<4 or len(d)<10:raise ValueError('Insufficient donors or assessable endpoints after explicit prerequisite exclusions')
    return d,spec


def dispersion_result(d,bootstrap=999,workers=1):
    def fit(sub):
        X=pd.DataFrame({'Intercept':1.,'disease':sub.donor_status.eq('T1D').astype(float),'age':sub.age,
                        'sex':sub.sex.astype(str).str.lower().str.startswith('m').astype(float)},index=sub.index)
        if np.linalg.matrix_rank(X)<X.shape[1]:raise ValueError('Rank-deficient donor dispersion design')
        return float(sm.OLS(sub.response,X).fit().params['disease'])
    value=fit(d);rng=np.random.default_rng(SEED)
    def one(ids):
        try:return fit(d.iloc[ids])
        except ValueError:return np.nan
    draws=np.array(Parallel(n_jobs=workers)(delayed(one)(rng.integers(0,len(d),len(d))) for _ in range(bootstrap)));draws=draws[np.isfinite(draws)]
    if len(draws)<.8*bootstrap:raise ValueError('Too few successful donor-level dispersion bootstrap fits')
    ci=np.quantile(draws,[.025,.975]);se=np.std(draws,ddof=1)
    return dict(effect=value,ci_low=float(ci[0]),ci_high=float(ci[1]),p_value=float(2*norm.sf(abs(value/se))) if se>0 else 1.,bootstrap_success=len(draws))


def prediction_result(base,regions,profiles,bootstrap,workers,root):
    d=base.loc[base.donor_status.isin(['ND','T1D'])].merge(region_wide(regions),on='islet_id',validate='one_to_one')
    if 'primary_assessable' in d:d=d.loc[d.primary_assessable.fillna(False)]
    d=d.merge(profiles.loc[profiles.profile_eligible],on=['case_id','islet_id'],validate='one_to_one')
    structural=['nerve_area_fraction_core','capillary_area_fraction_core','lymphatic_area_fraction_peri','coliv_boundary_coverage_core','duct_distance_um_core']
    if any(d[c].isna().all() for c in structural):raise ValueError('Full H10 profile lacks required validated structural/spatial measurements')
    coupling=weighted_coupling(pd.read_parquet(root/'analysis/coupling.parquet'),20)
    d=d.merge(coupling[['islet_id','excess_proximity']],on='islet_id',validate='one_to_one');structural+=['excess_proximity']
    d['male']=d.sex.astype(str).str.lower().str.startswith('m').astype(float)
    components=[c for c in d if c.startswith('endocrine_component__')]
    if len(components)!=4:raise ValueError('H10 requires all four expanded endocrine composition components')
    fractions=[]
    for column in components:
        name=column+'_fraction';d[name]=d[column]/d.trials;fractions.append(name)
    baseline=['age','male','log_area']+fractions;full=baseline+MARKERS+structural
    donor=d.groupby(['case_id','donor_status'])[full].mean().reset_index()
    if donor.case_id.nunique()<10 or donor[baseline].isna().any().any():raise ValueError('H10 baseline unavailable in sufficient donors')
    result=permutation_test(donor[baseline].to_numpy(float),donor[full].to_numpy(float),donor.donor_status.eq('T1D').to_numpy(int),bootstrap=bootstrap,workers=workers,exact=True)
    result['feature_names']=full
    return donor,result


def coefficient_stability(coefficients,features):
    coefficients=np.asarray(coefficients,float)
    if coefficients.ndim!=2 or coefficients.shape[1]!=len(features):raise ValueError('Feature stability coefficient dimensions differ')
    n=coefficients.shape[0];top=max(1,int(np.ceil(len(features)/4)))
    selected=np.zeros(coefficients.shape,bool)
    for i,row in enumerate(coefficients):
        order=sorted(range(len(features)),key=lambda j:(-abs(row[j]),features[j]))
        selected[i,order[:top]]=True
    rows=[]
    for j,name in enumerate(features):
        sign=max(float(np.mean(coefficients[:,j]>0)),float(np.mean(coefficients[:,j]<0))) if n else 0
        frequency=float(selected[:,j].mean()) if n else 0
        rows.append(dict(feature=name,sign_agreement=sign,top_quartile_frequency=frequency,
                         stable=bool(n>=10 and sign>=.8 and frequency>=.8 and name!='INS' and (name in MARKERS or name in ['nerve_area_fraction_core','capillary_area_fraction_core','lymphatic_area_fraction_peri','coliv_boundary_coverage_core','duct_distance_um_core'])),replicates=n))
    return rows


def dispersion_feature_stability(d,bootstrap=999):
    cols=['dispersion__'+m for m in MARKERS];rng=np.random.default_rng(SEED);coefficients=[]
    for _ in range(bootstrap):
        sample=d.iloc[rng.integers(0,len(d),len(d))]
        X=np.c_[np.ones(len(sample)),sample.donor_status.eq('T1D').astype(float),sample.age,sample.sex.astype(str).str.lower().str.startswith('m').astype(float)]
        Y=sample[cols].to_numpy(float);scale=Y.std(axis=0);scale[scale==0]=1
        if np.linalg.matrix_rank(X)<X.shape[1]:continue
        coefficients.append(np.linalg.lstsq(X,Y/scale,rcond=None)[0][1])
    return coefficient_stability(np.asarray(coefficients),MARKERS)


def followups(results,tables,specs,bootstrap=999,workers=1,root=None,base=None,regions=None,donors=None,report_root=None):
    r=results.set_index('id');estimated=set(r.index[r.p_value.notna()]);supported=set(r.index[r.q_value<.05]);rows=[]
    def missing(fid,required):return dict(id=fid,status='unavailable',reason='Required primary endpoints unavailable: '+', '.join(sorted(required-estimated)))
    if {'H1','H2','H3'}<=estimated:
        chosen=min(['H2','H3'],key=lambda h:(r.loc[h,'q_value'],0 if h=='H3' else 1))
        modify='H1' in supported and bool({'H2','H3'}&supported);d=tables['H1'].copy()
        d['nerve']=d.nerve_area_fraction_core;d['capillary']=d.capillary_area_fraction_core
        spec=dict(family='beta-binomial',term='disease:log_area_within',interaction=True,extra=['nerve','capillary'])
        if modify:spec.update(term='product',product=('log_area_within',('nerve' if chosen=='H2' else 'capillary')+'_within'))
        rec=dict(id='F1',branch=('size_by_'+chosen if modify else 'H1_adjust_both_structures'),status='exploratory_same_cohort')
        try:rec.update(bootstrap_model(d.dropna(subset=['nerve','capillary']),spec,bootstrap,workers))
        except (ValueError,np.linalg.LinAlgError) as e:rec.update(status='unavailable',reason=str(e))
        rows.append(rec)
    else:rows.append(missing('F1',{'H1','H2','H3'}))
    # These branches require explicitly exported neighborhoods/held-out
    # signatures; unavailable primary results never select a negative branch.
    for fid,required in [('F2',{'H4'}),('F3',{'H7','H8'}),('F4',{'H9','H10'})]:
        if not required<=estimated:rows.append(missing(fid,required));continue
        if fid=='F3':
            active=required&supported
            if not active:
                rows.append(dict(id=fid,status='unsupported_mechanism',branch='neither_primary_supported',reason='Neither adjusted primary test supports the prespecified mechanism; intervals retained in primary tables.'));continue
            h='H7' if 'H7' in active else 'H8';d=tables[h].copy();spec=dict(specs[h])
            if len(active)==2:
                d['duct_distance']=d.duct_distance_um_core;spec.update(extra=['duct_distance'],product=('coliv_within','duct_distance_within'),term='product');branch='duct_by_coliv'
            else:spec.update(disease_interaction=True,term='disease:exposure');branch=h+'_disease_interaction'
            try:rows.append(dict(id=fid,status='exploratory_same_cohort',branch=branch,**bootstrap_model(d,spec,bootstrap,workers)))
            except (ValueError,np.linalg.LinAlgError) as e:rows.append(dict(id=fid,status='unavailable',branch=branch,reason=str(e)))
        elif fid=='F2':
            from exploratory_v2 import lymphatic_subsets,lymphatic_localization
            branch='lymphatic_subsets' if 'H4' in supported else 'peri_vs_distant_lymphatics'
            rec=dict(id=fid,branch=branch,status='exploratory_same_cohort')
            try:
                if root is None:raise ValueError('Explicit endpoint root required for lymphatic follow-up')
                if 'H4' in supported:
                    data=pd.read_parquet(root/'analysis/lymphatic_neighborhoods.parquet')
                    if 'validated' not in data:raise ValueError('Neighborhood endpoints require affirmative validation')
                    data=data.loc[data.validated.fillna(False)].merge(base[['islet_id','case_id','age','sex','log_area','donor_status']],on=['islet_id','case_id'],validate='many_to_one')
                    rec.update(lymphatic_subsets(data,bootstrap,workers))
                else:rec.update(lymphatic_localization(pd.read_parquet(root/'analysis/lymphatic_localization.parquet'),bootstrap))
            except (ValueError,KeyError,FileNotFoundError,np.linalg.LinAlgError) as e:rec.update(status='unavailable',reason=str(e))
            rows.append(rec)
        else:
            from exploratory_v2 import residual_niche
            rec=dict(id=fid,branch='held_out_residual_endocrine_signature',status='exploratory_same_cohort')
            try:
                if root is None:raise ValueError('Explicit endpoint root required for residual niche follow-up')
                stable=root/'analysis/signature_stability.json'
                result,predictions=residual_niche(root,base,region_wide(regions),donors,bootstrap,workers,
                    stable_features=stable.exists() and json.loads(stable.read_text()).get('stable'),cache_root=report_root/'F4-reconstruction' if report_root else root/'analysis/F4-reconstruction')
                rec.update(result)
                if report_root:predictions.to_csv(report_root/'F4-held-out-predictions.csv',index=False)
            except (ValueError,KeyError,FileNotFoundError,np.linalg.LinAlgError) as e:rec.update(status='unavailable',reason=str(e))
            rows.append(rec)
    return pd.DataFrame(rows)


def write_plots(output,results,tables):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    good=results.loc[results.effect.notna()&results.ci_low.notna()&results.ci_high.notna()]
    if len(good):
        fig,ax=plt.subplots(figsize=(8,max(3,len(good)*.5)))
        for j,row in enumerate(good.itertuples()):ax.plot([row.ci_low,row.ci_high],[j,j],color='black');ax.plot(row.effect,j,'o')
        ax.set_yticks(np.arange(len(good)),good.id);ax.axvline(0,color='gray',lw=.8);ax.set_xlabel('Primary estimate (endpoint-specific units)');fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(output/f'primary-effects.{ext}',dpi=160)
        plt.close(fig)
    for h,d in tables.items():
        if 'islet_id' not in d or not len(d):continue
        response=d.success/d.trials if {'success','trials'}<=set(d) and h not in ['H2','H3','H4','H5'] else d.get('response')
        if response is None:continue
        fig,ax=plt.subplots(figsize=(7,4))
        for status,g in d.assign(display_response=response).groupby('donor_status'):ax.scatter(g.log_area,g.display_response,s=5,alpha=.2,label=status)
        ax.set(xlabel='log section area (µm²)',ylabel='Observed endpoint',title=h+' — descriptive observations');ax.legend();fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(output/f'{h}-observations.{ext}',dpi=150)
        plt.close(fig)
        d.assign(response_display=response,distance=(response-response.median()).abs()).sort_values(['distance','islet_id']).groupby('donor_status').head(1).to_csv(output/f'{h}-representative-islets.csv',index=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--report-root',required=True)
    p.add_argument('--registry',default=str(Path(__file__).with_name('hypotheses-v2.json')));p.add_argument('--bootstrap',type=int,default=999);p.add_argument('--workers',type=int,default=1);a=p.parse_args()
    root=Path(a.input_root);out=offline_report_root(root,a.report_root);registry=json.loads(Path(a.registry).read_text())
    if (out/'results.csv').exists():raise ValueError('Use a fresh report root; existing results are immutable')
    base,regions,summary,profiles,donors=prepare(root);profiles.to_parquet(root/'analysis/endocrine_profiles.parquet',index=False)
    memberships=pd.read_parquet(root/'analysis/memberships.parquet');memberships.case_id=memberships.case_id.astype(str)
    joint_counts(root,memberships,donors.case_id).to_parquet(root/'analysis/joint_marker_counts.parquet',index=False)
    results=[];tables={};specs={};sensitivities=[];stability=[]
    for h in registry['hypotheses']:
        hid=h['id'];rec=dict(id=hid,title=h['title'],status='unavailable',reason='',effect=np.nan,ci_low=np.nan,ci_high=np.nan,p_value=np.nan,n_donors=0,n_islets=0)
        try:
            if hid=='H10':
                d,result=prediction_result(base,regions,profiles,a.bootstrap,a.workers,root)
                pd.DataFrame({'null_delta':result['null']}).to_csv(out/'H10-permutations.csv',index=False)
                stability.extend([dict(hypothesis='H10',**record) for record in coefficient_stability(result['coefficients'],result['feature_names'])])
                rec.update({k:v for k,v in result.items() if np.isscalar(v)})
            else:
                d,spec=endpoint(hid,base,regions,summary,profiles,root);tables[hid]=d;specs[hid]=spec
                if hid=='H9':stability.extend([dict(hypothesis='H9',**record) for record in dispersion_feature_stability(d,a.bootstrap)])
                rec.update(dispersion_result(d,a.bootstrap,a.workers) if hid=='H9' else bootstrap_model(d,spec,a.bootstrap,a.workers,checkpoint_root=out/'checkpoints'/hid))
                if hid!='H9':
                    for label,band,flagged in [('band10',10,False),('band30',30,False),('include_flagged',20,True)]:
                        if hid=='H1' and band!=20:continue
                        try:
                            sd,ss=endpoint(hid,base,regions,summary,profiles,root,band,flagged)
                            result=bootstrap_model(sd,ss,a.bootstrap,a.workers,checkpoint_root=out/'checkpoints'/hid/label)
                            sensitivities.append(dict(id=hid,analysis=label,n_donors=sd.case_id.nunique(),n_islets=len(sd),**result))
                        except (ValueError,KeyError,FileNotFoundError,np.linalg.LinAlgError) as e:sensitivities.append(dict(id=hid,analysis=label,error=str(e)))
                    if 'donor_source' in d:
                        sd=d.loc[d.donor_source.ne('HPAP')]
                        try:sensitivities.append(dict(id=hid,analysis='exclude_HPAP',n_donors=sd.case_id.nunique(),**bootstrap_model(sd,spec,a.bootstrap,a.workers,checkpoint_root=out/'checkpoints'/hid/'exclude_HPAP')))
                        except (ValueError,np.linalg.LinAlgError) as e:sensitivities.append(dict(id=hid,analysis='exclude_HPAP',error=str(e)))
                for donor in sorted(d.case_id.unique()):
                    sub=d.loc[d.case_id!=donor]
                    try:value=dispersion_result(sub,max(20,min(99,a.bootstrap)),1)['effect'] if hid=='H9' else estimate(sub,**spec);sensitivities.append(dict(id=hid,analysis='omit_'+donor,effect=value,n_donors=sub.case_id.nunique()))
                    except ValueError as e:sensitivities.append(dict(id=hid,analysis='omit_'+donor,error=str(e)))
            rec.update(status='estimated_pending_scientific_review',n_donors=int(d.case_id.nunique()),n_islets=int(d.islet_id.nunique()) if 'islet_id' in d else 0)
            d.to_csv(out/f'{hid}-endpoint.csv',index=False)
        except (ValueError,KeyError,FileNotFoundError,np.linalg.LinAlgError) as e:rec['reason']=str(e)
        results.append(rec);print(hid,rec['status'],rec['reason'],flush=True)
    results=pd.DataFrame(results);results['q_value']=bh_ten(results.p_value);results['primary_family']='ten_hypotheses_v2';results.to_csv(out/'results.csv',index=False)
    write_json(root/'analysis/signature_stability.json',dict(stable=any(r['stable'] for r in stability),features=stability,
        definition='At least 80% same-sign and 80% top-quartile absolute standardized coefficients over >=10 reconstructed fits; INS excluded from signature candidates',
        marker_summaries_sha256=sha256(root/'analysis/marker_summaries.parquet')))
    followups(results,tables,specs,a.bootstrap,a.workers,root,base,regions,donors,out).to_csv(out/'followups.csv',index=False)
    pd.DataFrame(sensitivities).to_csv(out/'sensitivities.csv',index=False);write_plots(out,results,tables)
    write_json(out/'registry.json',registry)
    write_json(out/'run.json',dict(version=2,registry_sha256=sha256(a.registry),runner_sha256=sha256(__file__),
               inputs={str(p.relative_to(root)):sha256(p) for p in [root/'analysis/islets.parquet',root/'analysis/regions.parquet',root/'analysis/marker_summaries.parquet',root/'trajectory/endocrine_counts.parquet']},
               completed_tests=int(results.p_value.notna().sum()),bootstrap_requested=a.bootstrap,scientific_acceptance=False,
               acquisition='One acquisition per donor; acquisition and donor effects are inseparable.',pseudotime_used=False))

if __name__=='__main__':main()
