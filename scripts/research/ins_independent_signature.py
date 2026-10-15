"""Nested donor-held-out INS-independent signature with reconstructed cell typing.

INS is omitted before cell discovery, scaling and projection. Outer and inner
folds each get their own typing model and typed endocrine profiles. Conditional
training-reference identities never replace the validated production classifier.
"""
from __future__ import annotations
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel,delayed
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
from common import MARKERS,SEED,sha256,write_json
from cell_typing import TypingConfig,discover,read_cells,project_frame,consensus_assignments,atomic_parquet,canonical_hash,attach_quality_status
from cell_typing_validation import align_training_clusters
from conditional_typing import hierarchy_reference,reconstruct as reconstruct_typing_hierarchy
from exploratory_v2 import interval

FEATURE_MARKERS=[m for m in MARKERS if m!='INS']
STRUCTURAL=['log_area','nerve_area_fraction_core','capillary_area_fraction_core','lymphatic_area_fraction_peri']


def reconstruct_profiles(root,reference,cache,excluded,workers=1,fit_workers=1):
    excluded=sorted(set(excluded));_,_,source,hierarchy_identity=hierarchy_reference(reference,excluded)
    folder=cache/('without_'+'_'.join(excluded));checkpoint=folder/'profile_complete.json'
    donors=pd.read_csv(root/'cohort.csv',dtype=str).case_id.tolist()
    fingerprint=canonical_hash(dict(excluded=excluded,corrected=source['training_inputs'],typing_config=source['config'],
        labels={c:sha256(reference/'typing/cells'/f'{c}.parquet') for c in donors if c not in excluded},
        hierarchy=hierarchy_identity,membership=sha256(root/'analysis/memberships.parquet'),code=sha256(__file__)))
    if checkpoint.exists():
        if json.loads(checkpoint.read_text())['fingerprint']!=fingerprint:raise ValueError('INS-independent profile cache changed inputs')
        return pd.read_parquet(folder/'profiles.parquet')
    spec,config,predictions=reconstruct_typing_hierarchy(root,folder,reference,excluded,{'INS'},donors,workers=workers,fit_workers=fit_workers,endocrine_only=True)
    members=pd.read_parquet(root/'analysis/memberships.parquet');members.case_id=members.case_id.astype(str);members.object_id=members.object_id.astype(str)
    members=members.loc[members.compartment.eq('core')&members.band_um.eq(20)]
    profiles=[]
    for cid in donors:
        cells=read_cells(root/'corrected'/f'{cid}-cells.parquet',cid)
        labels=pd.read_parquet(predictions/f'{cid}.parquet')
        lookup=pd.Index(labels.object_id).get_indexer(cells.object_id)
        if (lookup<0).any() or len(lookup)!=len(labels):raise ValueError('Conditional endocrine UUID universe differs from raw cells')
        labels=labels.iloc[lookup].reset_index(drop=True)
        endocrine=labels.is_endocrine.fillna(False)&labels.classification_status.eq('diagnostic_label_transfer')&labels.projection_supported&~labels.segmentation_quality_excluded
        selected=cells.loc[endocrine,['case_id','object_id']+FEATURE_MARKERS].copy()
        selected[FEATURE_MARKERS]=np.log1p(selected[FEATURE_MARKERS].to_numpy(float))
        joined=members.loc[members.case_id.eq(cid)].merge(selected,on=['case_id','object_id'],validate='many_to_one')
        grouped=joined.groupby(['islet_id','case_id']);means=grouped[FEATURE_MARKERS].mean();means['endocrine_measured_cells']=grouped.size()
        profiles.append(means.loc[means.endocrine_measured_cells>=5].reset_index())
    frame=pd.concat(profiles,ignore_index=True);atomic_parquet(frame,folder/'profiles.parquet')
    write_json(checkpoint,dict(fingerprint=fingerprint,excluded_donors=excluded,omitted_markers=['INS'],
        training_donors=sorted(spec['training_inputs']),typing_status='conditional_training_reference_projection',
        minimum_endocrine_cells=5,features=FEATURE_MARKERS,
        interpretation='INS-independent predictor measurements; broader endocrine-module pseudotime exclusion is a separate analysis.'))
    return frame


def select_features(training,columns,maximum=8):
    # Candidate support is fixed, but selection uses training donors only.
    donor=training.groupby('case_id')[columns+['insulin_burden']].mean()
    score=[]
    for column in columns:
        v=donor[[column,'insulin_burden']].dropna()
        r=spearmanr(v[column],v.insulin_burden).statistic if len(v)>=4 and v[column].nunique()>1 else np.nan
        if np.isfinite(r):score.append((abs(float(r)),column))
    if not score:raise ValueError('No supported variable features in INS-independent training donors')
    return [c for _,c in sorted(score,key=lambda pair:(-pair[0],pair[1]))[:maximum]]


def resample_training(frame,draw):
    if draw is None:return frame
    parts=[]
    for j,cid in enumerate(draw):
        g=frame.loc[frame.case_id.eq(cid)].copy();g['source_case_id']=g.case_id;g.case_id=f'{cid}@{j}';parts.append(g)
    return pd.concat(parts,ignore_index=True) if parts else frame.iloc[:0]


def fit_ridge(train,test,columns,alpha):
    weight=1/train.groupby('case_id').case_id.transform('size').to_numpy()
    model=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),Ridge(alpha=alpha))
    model.fit(train[columns],train.insulin_burden,ridge__sample_weight=weight)
    return model.predict(test[columns])


def predict_fold(held,frames,endpoint,donors,draw=None):
    population=sorted(set(draw if draw is not None else donors))
    train_donors=[c for c in population if c!=held]
    if len(train_donors)<4:raise ValueError('Insufficient distinct training donors for nested signature')
    columns=FEATURE_MARKERS+STRUCTURAL
    def data(excluded):
        key=tuple(sorted(excluded));f=frames[key].merge(endpoint,on=['islet_id','case_id'],validate='one_to_one')
        return f.dropna(subset=['insulin_burden']+STRUCTURAL)
    outer=data([held]);train=outer.loc[outer.case_id.isin(train_donors)];test=outer.loc[outer.case_id.eq(held)]
    train=resample_training(train,None if draw is None else [c for c in draw if c!=held])
    if len(test)<3:raise ValueError('Insufficient supported endocrine profiles in held-out donor')
    candidates=[.01,.1,1.,10.];scores=[]
    for alpha in candidates:
        errors=[]
        for inner in train_donors:
            frame=data([held,inner]);validation=frame.loc[frame.case_id.eq(inner)]
            fit=frame.loc[frame.case_id.isin([c for c in train_donors if c!=inner])]
            fit=resample_training(fit,None if draw is None else [c for c in draw if c not in {held,inner}])
            if not len(validation) or fit.case_id.nunique()<3:continue
            chosen=select_features(fit,columns);prediction=fit_ridge(fit,validation,chosen,alpha)
            errors.append(float(np.mean((prediction-validation.insulin_burden.to_numpy())**2)))
        if len(errors)<.8*len(train_donors):raise ValueError('Too few supported inner folds for INS-independent signature')
        scores.append(float(np.mean(errors)))
    alpha=candidates[int(np.argmin(scores))];chosen=select_features(train,columns)
    result=test[['islet_id','case_id','insulin_burden']].copy();result['held_out_signature']=fit_ridge(train,test,chosen,alpha)
    return result,dict(held_out=held,alpha=alpha,features=chosen,inner_mse=scores,supported_islets=len(result))


def run_signature(root,reference,cache,endpoint,workers=1,fit_workers=1,bootstrap=999):
    donors=sorted(endpoint.case_id.unique());frames={};predictions=[];audits=[]
    # Pairwise excluded sets are shared by both outer/inner orientations.
    needed={(d,) for d in donors}|{tuple(sorted((a,b))) for i,a in enumerate(donors) for b in donors[i+1:]}
    for excluded in sorted(needed,key=lambda x:(len(x),x)):
        frames[excluded]=reconstruct_profiles(root,reference,cache/'typing',excluded,workers,fit_workers)
    for held in donors:
        prediction,audit=predict_fold(held,frames,endpoint,donors);predictions.append(prediction);audits.append(audit)
    predictions=pd.concat(predictions,ignore_index=True)
    means=predictions.groupby('case_id')[['held_out_signature','insulin_burden']].mean()
    effect=float(spearmanr(means.held_out_signature,means.insulin_burden).statistic)
    rng=np.random.default_rng(SEED);draws=[rng.choice(donors,len(donors),replace=True).tolist() for _ in range(bootstrap)]
    checkpoint=cache/'bootstrap';checkpoint.mkdir(parents=True,exist_ok=True)
    fingerprint=hashlib.sha256(pd.util.hash_pandas_object(endpoint,index=False).to_numpy().tobytes()+sha256(__file__).encode()+b''.join(sha256(cache/'typing'/('without_'+'_'.join(k))/'profiles.parquet').encode() for k in sorted(frames))).hexdigest()
    def one(index,draw):
        path=checkpoint/f'{index}.json'
        if path.exists():
            saved=json.loads(path.read_text())
            if saved['fingerprint']==fingerprint:return np.nan if saved['value'] is None else saved['value']
        pairs=[]
        for held in sorted(set(draw)):
            try:
                with threadpool_limits(limits=1):p,_=predict_fold(held,frames,endpoint,donors,draw)
                pair=[float(p.held_out_signature.mean()),float(p.insulin_burden.mean())];pairs.extend([pair]*draw.count(held))
            except (ValueError,np.linalg.LinAlgError):pass
        value=float(spearmanr(np.asarray(pairs)[:,0],np.asarray(pairs)[:,1]).statistic) if len(pairs)>=.8*len(draw) else np.nan
        write_json(path,dict(fingerprint=fingerprint,value=value if np.isfinite(value) else None));return value
    values=Parallel(n_jobs=workers)(delayed(one)(i,draw) for i,draw in enumerate(draws))
    result=dict(interval(effect,values,bootstrap),branch='stable_INS_independent_nested_signature',
        unit='Donor-mean correlation of held-out INS-independent signature with insulin-positive tissue burden',
        typing_omitted_markers=['INS'],n_donors=len(donors),n_islets=len(predictions),feature_selection='Top eight absolute training-donor Spearman associations; repeated inside inner and outer folds',
        uncertainty='Donor bootstrap repeats nested feature selection/regression on cached independently reconstructed fold typing; production identity uncertainty remains conditional on training references.')
    write_json(cache/'fold_models.json',audits);atomic_parquet(predictions,cache/'held_out_predictions.parquet')
    return result,predictions
