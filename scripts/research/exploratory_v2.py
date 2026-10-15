"""Result-dependent exploratory follow-ups with explicit endpoint contracts.

These tests are not additional primary hypotheses. Data extraction and marker
review must finish first; missing inputs never choose an unsupported mechanism.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2,norm,spearmanr
from scipy.optimize import minimize
from joblib import Parallel,delayed
from threadpoolctl import threadpool_limits
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from common import SEED
from run_hypotheses import beta_binomial_fit,design,within_between


def donor_resamples(data,repetitions,seed=SEED):
    donors=np.array(sorted(data.case_id.unique()));rng=np.random.default_rng(seed)
    for _ in range(repetitions):
        parts=[]
        for j,cid in enumerate(rng.choice(donors,len(donors),replace=True)):
            part=data.loc[data.case_id.eq(cid)].copy();part['source_case_id']=cid;part.case_id=f'{cid}@{j}';parts.append(part)
        yield pd.concat(parts,ignore_index=True)


def interval(effect,draws,requested):
    values=np.asarray(draws,float);values=values[np.isfinite(values)]
    if len(values)<max(2,.8*requested):raise ValueError('Fewer than 80% exploratory donor bootstrap fits succeeded')
    ci=np.quantile(values,[.025,.975]);se=np.std(values,ddof=1)
    return dict(effect=float(effect),ci_low=float(ci[0]),ci_high=float(ci[1]),p_value=float(2*norm.sf(abs(effect/se))) if se>0 else 1.,bootstrap_success=len(values))


def lymphatic_subsets(data,bootstrap=999,workers=1):
    subsets=['CD8','CD4','B','macrophage']
    d=data.loc[data.total_count.gt(0)].dropna(subset=['near_count','total_count','near_area_fraction','age','sex','log_area']).copy()
    if set(d.subset)!=set(subsets):raise ValueError('All four validated immune subset neighborhood endpoints are required')
    def fit(sub):
        X=design(sub)
        coverage=within_between(sub,'near_area_fraction','coverage');X['coverage_within']=coverage.coverage_within;X['coverage_between']=coverage.coverage_between
        terms=[]
        for subset in subsets[1:]:
            term='subset_'+subset;X[term]=sub.subset.eq(subset).astype(float);terms.append(term)
        if np.linalg.matrix_rank(X)<X.shape[1]:raise ValueError('Rank-deficient lymphatic subset model')
        weights=1/sub.groupby('case_id').case_id.transform('size').to_numpy();weights*=len(weights)/weights.sum()
        beta=beta_binomial_fit(X,sub.near_count,sub.total_count,weights)
        return beta[[list(X).index(t) for t in terms]]
    observed=fit(d)
    def one(sub):
        try:
            with threadpool_limits(limits=1):return fit(sub)
        except (ValueError,np.linalg.LinAlgError):return np.full(3,np.nan)
    draws=np.asarray(Parallel(n_jobs=workers)(delayed(one)(sub) for sub in donor_resamples(d,bootstrap)))
    draws=draws[np.isfinite(draws).all(axis=1)]
    if len(draws)<.8*bootstrap:raise ValueError('Too few successful subset omnibus bootstrap fits')
    covariance=np.cov(draws,rowvar=False)
    if np.linalg.matrix_rank(covariance)<3:raise ValueError('Singular donor bootstrap subset covariance')
    statistic=float(observed@np.linalg.solve(covariance,observed))
    return dict(effect=statistic,statistic='donor-bootstrap covariance Wald chi-square',df=3,p_value=float(chi2.sf(statistic,3)),
                contrasts=['CD4 versus CD8','B versus CD8','macrophage versus CD8'],log_odds=observed.tolist(),
                ci_low=np.quantile(draws,.025,axis=0).tolist(),ci_high=np.quantile(draws,.975,axis=0).tolist(),
                bootstrap_success=len(draws),n_donors=d.case_id.nunique(),n_islets=d.islet_id.nunique())


def lymphatic_localization(data,bootstrap=999):
    required={'case_id','compartment','positive_area_um2','valid_area_um2','validated'}
    if not required.issubset(data):raise ValueError('Missing peri/distant lymphatic localization measurements')
    d=data.loc[data.validated.fillna(False)&data.valid_area_um2.gt(0)].copy()
    if not d.compartment.isin(['peri','distant']).all():raise ValueError('Localization requires disjoint peri and distant tissue regions')
    grouped=d.groupby(['case_id','compartment'])[['positive_area_um2','valid_area_um2']].sum()
    grouped['fraction']=grouped.positive_area_um2/grouped.valid_area_um2
    paired=grouped.fraction.unstack().dropna();delta=(paired.peri-paired.distant).to_numpy(float)
    if len(delta)<4:raise ValueError('Insufficient donors with both assessed localization compartments')
    rng=np.random.default_rng(SEED);draws=[float(np.mean(rng.choice(delta,len(delta),replace=True))) for _ in range(bootstrap)]
    return dict(interval(float(delta.mean()),draws,bootstrap),n_donors=len(delta),unit='peri minus distant lymphatic positive-area fraction')


def nested_signature(data,columns):
    """One prediction for each donor-held-out islet, train-only preprocessing."""
    d=data.copy();groups=d.get('source_case_id',d.case_id).astype(str).to_numpy();X=d[columns].to_numpy(float);y=d.insulin_burden.to_numpy(float)
    prediction=np.full(len(d),np.nan)
    def model(alpha):return make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),Ridge(alpha=alpha))
    def weights(indices):
        _,inverse,n=np.unique(groups[indices],return_inverse=True,return_counts=True)
        return 1/n[inverse]
    for held in sorted(set(groups)):
        train=np.flatnonzero(groups!=held);test=np.flatnonzero(groups==held)
        if len(set(groups[train]))<3:raise ValueError('Too few donors for held-out signature reconstruction')
        errors=[]
        for alpha in [.01,.1,1.,10.]:
            donor_errors=[]
            for inner in sorted(set(groups[train])):
                fit_indices=train[groups[train]!=inner];valid=train[groups[train]==inner]
                fitted=model(alpha).fit(X[fit_indices],y[fit_indices],ridge__sample_weight=weights(fit_indices))
                donor_errors.append(float(np.mean((fitted.predict(X[valid])-y[valid])**2)))
            errors.append(np.mean(donor_errors))
        fitted=model([.01,.1,1.,10.][int(np.argmin(errors))]).fit(X[train],y[train],ridge__sample_weight=weights(train))
        prediction[test]=fitted.predict(X[test])
    d['held_out_signature']=prediction
    # Donors are biological replicates for the supporting association.
    donor=d.groupby('case_id')[['held_out_signature','insulin_burden']].mean()
    r=float(spearmanr(donor.held_out_signature,donor.insulin_burden).statistic)
    return r,d


def censored_c_peptide(data,bootstrap=999,workers=1):
    """Left-censored log-normal donor model; '<' direction is the metadata contract."""
    d=data.dropna(subset=['C_peptide_ng_mL','held_out_signature','age','sex']).copy()
    d=d.loc[d.C_peptide_ng_mL>0]
    if len(d)<8:raise ValueError('Insufficient donors for secondary censored C-peptide regression')
    def fit(sub):
        X=np.c_[np.ones(len(sub)),sub.held_out_signature,sub.age,
                sub.sex.astype(str).str.lower().str.startswith('m').astype(float),
                sub.donor_status.eq('T1D').astype(float),sub.donor_status.eq('Aab+').astype(float)]
        if np.linalg.matrix_rank(X)<X.shape[1]:raise ValueError('Rank-deficient censored C-peptide design')
        center=X.mean(axis=0);scale=X.std(axis=0);center[scale==0]=0;scale[scale==0]=1;Z=(X-center)/scale
        y=np.log(sub.C_peptide_ng_mL.to_numpy(float));censored=sub.C_peptide_censored.astype(str).str.lower().eq('true').to_numpy()
        def objective(theta):
            sigma=np.exp(theta[-1]);z=(y-Z@theta[:-1])/sigma
            ll=np.where(censored,norm.logcdf(z),norm.logpdf(z)-np.log(sigma))
            return -float(np.sum(ll))
        start=np.r_[np.linalg.lstsq(Z,y,rcond=None)[0],0.]
        result=minimize(objective,start,method='L-BFGS-B',bounds=[(None,None)]*X.shape[1]+[(-8,8)])
        if not result.success or not np.isfinite(result.fun):raise ValueError('Censored C-peptide model did not converge')
        return float(result.x[1]/scale[1])
    observed=fit(d)
    def one(sub):
        try:return fit(sub)
        except (ValueError,np.linalg.LinAlgError):return np.nan
    draws=Parallel(n_jobs=workers)(delayed(one)(sub) for sub in donor_resamples(d,bootstrap))
    return dict(interval(observed,draws,bootstrap),n_donors=len(d),n_left_censored=int(d.C_peptide_censored.astype(str).str.lower().eq('true').sum()),
                unit='log C-peptide per unit held-out signature',censoring='Left censoring at the recorded positive numeric bound, per build_donor_metadata.parse_c_peptide',
                adjustment='age, sex, separate ND/Aab+/T1D categories',uncertainty='Donor bootstrap conditional on held-out signature predictions')


def residual_niche(root,base,region_table,donors,bootstrap=999,workers=1,stable_features=False,cache_root=None):
    path=root/'analysis/insulin_tissue_burden.parquet'
    if not path.exists():raise ValueError('Validated insulin-positive tissue area endpoint is missing; cell fractions cannot substitute for tissue burden')
    burden=pd.read_parquet(path)
    required={'islet_id','ins_positive_area_um2','valid_core_area_um2','validated'}
    if not required.issubset(burden):raise ValueError('Insulin tissue burden lacks its validated area denominator')
    burden=burden.loc[burden.validated.fillna(False)&burden.valid_core_area_um2.gt(0)].copy()
    burden['insulin_burden']=burden.ins_positive_area_um2/burden.valid_core_area_um2
    if not burden.insulin_burden.between(0,1).all():raise ValueError('Insulin tissue burden outside [0,1]')
    d=base.merge(region_table,on='islet_id',validate='one_to_one').merge(burden[['islet_id','insulin_burden']],on='islet_id',validate='one_to_one')
    # The fixed fallback is independent of all cell-type marker inputs. Stable
    # INS-excluded typed signatures require separately reconstructed held-out
    # artifacts; merely dropping INS from full-cohort typing is prohibited.
    columns=['log_area','nerve_area_fraction_core','capillary_area_fraction_core','lymphatic_area_fraction_peri']
    d=d.dropna(subset=columns+['insulin_burden'])
    if d.case_id.nunique()<5:raise ValueError('Insufficient fully assessed donors for the prespecified structural signature')
    if stable_features:
        from ins_independent_signature import run_signature
        cache=Path(cache_root) if cache_root is not None else root/'analysis/INS-independent-reconstruction'
        result,predictions=run_signature(root,root,cache,d,workers=workers,fit_workers=min(4,workers),bootstrap=bootstrap)
    else:
        result,predictions=fallback_niche_signature(d,columns,bootstrap,workers)
    clinical=predictions.groupby('case_id').held_out_signature.mean().reset_index().merge(donors,on='case_id',validate='one_to_one')
    try:
        result['c_peptide_status']='estimated_secondary_left_censored'
        result['c_peptide_association']=censored_c_peptide(clinical,bootstrap,workers)
    except ValueError as error:
        result['c_peptide_status']='unavailable';result['c_peptide_reason']=str(error)
    return result,predictions


def fallback_niche_signature(d,columns,bootstrap,workers):
    observed,predictions=nested_signature(d,columns)
    def one(sub):
        try:
            with threadpool_limits(limits=1):return nested_signature(sub,columns)[0]
        except (ValueError,np.linalg.LinAlgError):return np.nan
    draws=Parallel(n_jobs=workers)(delayed(one)(sub) for sub in donor_resamples(d,bootstrap))
    result=dict(interval(observed,draws,bootstrap),branch='prespecified_size_plus_structural_density',
                unit='Spearman correlation of donor-mean held-out signature with insulin-positive tissue burden',
                n_donors=d.case_id.nunique(),n_islets=d.islet_id.nunique(),features=columns)
    return result,predictions
