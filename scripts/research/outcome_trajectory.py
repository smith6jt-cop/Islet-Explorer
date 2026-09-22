#!/usr/bin/env python3
"""Outcome-module-excluded, donor-held-out supporting trajectory associations.

Each outer fold reconstructs all four typing fits without the outcome module,
using only training donors for scaling, discovery and cluster label alignment.
Labels transferred from reviewed training references are conditional analysis
labels; they never replace production cell identities. The trajectory is then
fit on the training islets and projected into the held-out donor. A smooth
association is evaluated by its held-out mean squared error improvement over
the same demographic/size covariates. No cross-fold axis alignment uses outcomes.
"""
from __future__ import annotations
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel,delayed
from sklearn.preprocessing import SplineTransformer
from scipy.stats import norm
from threadpoolctl import threadpool_limits
from common import MARKERS,SEED,sha256,write_json,offline_report_root
from composition import canonical_cells,aggregate_donor
from trajectory import MODULES,fit,project
from cell_typing import TypingConfig,discover,read_cells,project_frame,consensus_assignments,atomic_parquet,safe_name,configure_runtime,attach_quality_status
from cell_typing_validation import align_training_clusters
from conditional_typing import hierarchy_reference,reconstruct as reconstruct_typing_hierarchy

MODULE_LINEAGES=[{'endocrine'},{'neural'},{'vascular','smooth muscle'},{'lymphatic'},{'immune'},
                 {'stromal','epithelial'},{'technical'}]


def outcome_module(outcome,catalog):
    seed={outcome} if outcome in MARKERS else set()
    match=[entry for entry in catalog.get('base_types',[]) if entry['name']==outcome]
    if match:
        seed.update(match[0].get('marker_dependencies',[]))
        lineage=str(match[0].get('lineage','')).lower()
        for markers,lineages in zip(MODULES,MODULE_LINEAGES):
            if lineage in lineages:seed.update(markers)
    if not seed:raise ValueError('Outcome requires an explicit marker or audited phenotype dependency')
    omit=set(seed)
    for module in MODULES:
        if omit&module:omit.update(module)
    lineages=set().union(*(ls for ms,ls in zip(MODULES,MODULE_LINEAGES) if ms&omit))
    excluded_types={entry['name'] for entry in catalog['base_types'] if set(entry.get('marker_dependencies',[]))&omit or str(entry.get('lineage','')).lower() in lineages}
    if outcome not in MARKERS:excluded_types.add(outcome)
    return omit,excluded_types


def response_values(root,outcome,compartment='both'):
    if outcome in MARKERS:
        s=pd.read_parquet(root/'analysis/marker_summaries.parquet')
        s=s.loc[s.marker.eq(outcome)&s.population.eq('All cells')].copy()
        s=s.rename(columns={'positive_count':'count','positive_available':'available'})
    else:
        s=pd.read_parquet(root/'trajectory/counts.parquet');s=s.loc[s.phenotype.eq(outcome)].copy()
    selected=['core','peri'] if compartment=='both' else [compartment]
    s=s.loc[s.compartment.isin(selected)]
    rows=[]
    for iid,g in s.groupby('islet_id'):
        valid=len(g)==len(selected) and g.available.fillna(False).all() and g.endocrine_count.notna().all() and g.endocrine_count.gt(0).all()
        rows.append(dict(islet_id=iid,response=float(np.log1p(g['count'].sum()/g.endocrine_count.iloc[0])) if valid else np.nan))
    return pd.DataFrame(rows)


def reconstruct_fold(root,output,reference,held,omitted,excluded_types,workers=1,fit_workers=1):
    directory=output/'folds'/safe_name(held);checkpoint=directory/'composition_complete.json'
    _,_,source_spec,hierarchy_identity=hierarchy_reference(reference,[held])
    catalog=json.loads((reference/'typing/catalog.json').read_text())
    donors=pd.read_csv(root/'cohort.csv',dtype=str).case_id.tolist()
    fingerprint=hashlib.sha256(json.dumps(dict(held=held,omitted=sorted(omitted),excluded_types=sorted(excluded_types),
        typing_inputs=source_spec['training_inputs'],typing_config=source_spec['config'],reference_typing={c:sha256(reference/'typing/cells'/f'{c}.parquet') for c in donors if c!=held},
        hierarchy=hierarchy_identity,geometry=sha256(root/'analysis/islets.parquet'),memberships=sha256(root/'analysis/memberships.parquet'),code=sha256(__file__)),sort_keys=True).encode()).hexdigest()
    if checkpoint.exists():
        saved=json.loads(checkpoint.read_text())
        if saved['fingerprint']!=fingerprint:raise ValueError('Outcome fold inputs changed; use a fresh output root')
        return pd.read_parquet(directory/'composition.parquet')
    spec,config,predictions=reconstruct_typing_hierarchy(root,directory,reference,[held],omitted,donors,excluded_types,workers,fit_workers)
    base_types=[entry['name'] for entry in catalog['base_types'] if entry['name'] not in excluded_types]
    islets=pd.read_parquet(root/'analysis/islets.parquet');islets.case_id=islets.case_id.astype(str)
    members=pd.read_parquet(root/'analysis/memberships.parquet');members.case_id=members.case_id.astype(str)
    frames=[]
    for cid in donors:
        corrected=read_cells(root/'corrected'/f'{cid}-cells.parquet',cid)
        labels=pd.read_parquet(predictions/f'{cid}.parquet')
        omitted_identity=labels.base_type.isin(excluded_types)
        labels.loc[omitted_identity,['base_type','is_endocrine']]=None
        labels.loc[omitted_identity,'classification_status']='outcome_dependent_identity_excluded'
        labels.loc[labels.classification_status.eq('diagnostic_label_transfer'),'classification_status']='outcome_excluded_projection'
        atomic_parquet(labels,directory/'conditional_typing'/f'{cid}.parquet')
        cells=canonical_cells(labels,corrected,accepted_statuses=('outcome_excluded_projection',))
        features,_,_=aggregate_donor(islets.loc[islets.case_id.eq(cid)],members.loc[members.case_id.eq(cid)],cells,base_types,include_marker_summaries=False)
        features['typing_analysis_status']='outcome_excluded_projection';frames.append(features)
    frame=pd.concat(frames,ignore_index=True);atomic_parquet(frame,directory/'composition.parquet')
    write_json(checkpoint,dict(fingerprint=fingerprint,held_out_donor=held,omitted_markers=sorted(omitted),excluded_types=sorted(excluded_types),
        training_donors=sorted(spec['training_inputs']),model_identity=spec['identity'],
        typing_scope='Conditional training-reference transfer; no production identity changes',
        endocrine_denominator='Reconstructed from this outcome-excluded fold only; no original-denominator substitution'))
    return frame


def clinical_design(frame):
    return np.c_[np.ones(len(frame)),pd.to_numeric(frame.age),frame.sex.astype(str).str.lower().str.startswith('m').astype(float),
                 frame.donor_status.eq('T1D').astype(float),frame.donor_status.eq('Aab+').astype(float),np.log(frame.area_um2)]


def held_out_smooth(frame,held,response,donors,bootstrap_draw=None):
    d=frame.merge(response,on='islet_id',how='left',validate='one_to_one')
    metadata=donors[[c for c in ['case_id','age','sex','donor_status'] if c=='case_id' or c not in d]]
    d=d.merge(metadata,on='case_id',how='left',validate='many_to_one')
    train=d.loc[d.case_id.ne(held)].copy();test=d.loc[d.case_id.eq(held)].copy()
    if bootstrap_draw is not None:
        pieces=[]
        for j,cid in enumerate(bootstrap_draw):
            if cid==held:continue
            part=train.loc[train.case_id.eq(cid)].copy();part.case_id=f'{cid}@{j}';part.islet_id=part.islet_id+f'@{j}';pieces.append(part)
        train=pd.concat(pieces,ignore_index=True) if pieces else train.iloc[:0]
    train_coordinates,model=fit(train);test_coordinates=project(test,model)
    # Inference uses the diffusion coordinate itself.  The fixed reference rank
    # is a display/binning coordinate and may spread exact raw ties by an
    # identifier-only hash; that arbitrary within-tie order must never become
    # an inferential covariate.
    coordinate='pseudotime_raw'
    coordinate_columns=['islet_id','pseudotime_raw','pseudotime_percentile']
    train=train.merge(train_coordinates[coordinate_columns],on='islet_id',validate='one_to_one').dropna(subset=[coordinate,'response','age','sex'])
    test=test.merge(test_coordinates[coordinate_columns],on='islet_id',validate='one_to_one').dropna(subset=[coordinate,'response','age','sex'])
    if len(test)<4 or train.case_id.nunique()<3:raise ValueError('Insufficient supported observations in held-out smooth association')
    base_train=clinical_design(train);base_test=clinical_design(test)
    center=base_train.mean(axis=0);scale=base_train.std(axis=0);center[scale==0]=0;scale[scale==0]=1
    base_train=(base_train-center)/scale;base_test=(base_test-center)/scale
    spline=SplineTransformer(n_knots=4,degree=3,knots='quantile',include_bias=False,extrapolation='linear')
    smooth_train=spline.fit_transform(train[[coordinate]]);smooth_test=spline.transform(test[[coordinate]])
    weight=1/train.groupby('case_id').case_id.transform('size').to_numpy();weight=np.sqrt(weight)
    def prediction(X,Z):return Z@np.linalg.lstsq(X*weight[:,None],train.response.to_numpy()*weight,rcond=None)[0]
    null=prediction(base_train,base_test);full=prediction(np.c_[base_train,smooth_train],np.c_[base_test,smooth_test])
    error_null=(test.response.to_numpy()-null)**2;error_full=(test.response.to_numpy()-full)**2
    result=test[['islet_id','case_id','pseudotime_raw','pseudotime_percentile','response']].copy()
    result['baseline_prediction']=null;result['smooth_prediction']=full
    return float(np.mean(error_null-error_full)),result


def run(root,output,reference,outcome,workers=1,fit_workers=1,bootstrap=200,compartment='both'):
    catalog=json.loads((reference/'typing/catalog.json').read_text());omitted,excluded_types=outcome_module(outcome,catalog)
    if omitted&{'INS','GCG','SST','CHGA'}:
        result=dict(outcome=outcome,status='descriptive_only',reason='Endocrine module exclusion removes the evidence required to reconstruct the agreed endocrine denominator')
        write_json(output/'result.json',result);return result
    response=response_values(root,outcome,compartment)
    if response.response.notna().sum()<10:raise ValueError('Validated outcome measurements unavailable')
    donors=pd.read_csv(root/'donor_metadata.csv',dtype={'case_id':str});frames={};folds=[];predictions=[]
    module_key=hashlib.sha256(json.dumps([sorted(omitted),sorted(excluded_types)]).encode()).hexdigest()[:16]
    module_output=output.parent/'module_models'/module_key
    for held in sorted(donors.case_id):
        frame=reconstruct_fold(root,module_output,reference,held,omitted,excluded_types,workers,fit_workers);frames[held]=frame
        try:
            value,prediction=held_out_smooth(frame,held,response,donors);predictions.append(prediction)
            folds.append(dict(case_id=held,improvement=value,supported=len(prediction),status='supported'))
        except (ValueError,np.linalg.LinAlgError) as error:folds.append(dict(case_id=held,improvement=np.nan,supported=0,status=str(error)))
    fold_table=pd.DataFrame(folds);atomic_parquet(fold_table,output/'fold_results.parquet')
    if predictions:atomic_parquet(pd.concat(predictions,ignore_index=True),output/'held_out_predictions.parquet')
    supported=fold_table.improvement.notna().sum()
    if supported<.8*len(donors):raise ValueError('Fewer than 80% of donor folds support the reconstructed association')
    effect=float(fold_table.improvement.mean());rng=np.random.default_rng(SEED);draws=[rng.choice(sorted(frames),len(frames),replace=True) for _ in range(bootstrap)]
    checkpoints=output/'bootstrap';checkpoints.mkdir(parents=True,exist_ok=True)
    fingerprint=hashlib.sha256(pd.util.hash_pandas_object(response,index=False).to_numpy().tobytes()+sha256(__file__).encode()+b''.join(sha256(module_output/'folds'/safe_name(h)/'composition.parquet').encode() for h in sorted(frames))).hexdigest()
    def one(i,draw):
        path=checkpoints/f'{i}.json'
        if path.exists():
            record=json.loads(path.read_text())
            if record['fingerprint']==fingerprint:return np.nan if record['value'] is None else record['value']
        values=[]
        for held in sorted(set(draw)):
            try:
                with threadpool_limits(limits=1):value,_=held_out_smooth(frames[held],held,response,donors,draw)
                values.extend([value]*int(np.sum(draw==held)))
            except (ValueError,np.linalg.LinAlgError):pass
        value=float(np.mean(values)) if len(values)>=.8*len(draw) else np.nan
        write_json(path,dict(fingerprint=fingerprint,value=value if np.isfinite(value) else None));return value
    values=np.asarray(Parallel(n_jobs=workers)(delayed(one)(i,draw) for i,draw in enumerate(draws)));valid=values[np.isfinite(values)]
    complete=len(valid)>=.8*bootstrap;ci=np.quantile(valid,[.025,.975]) if len(valid)>=2 else [np.nan,np.nan];se=float(np.std(valid,ddof=1)) if len(valid)>=2 else np.nan
    result=dict(outcome=outcome,compartment=compartment,status='estimated_supporting_conditional_typing' if complete else 'incomplete_bootstrap',
        effect=effect,unit='Donor-mean held-out MSE improvement from cubic trajectory smooth',ci_low=float(ci[0]) if np.isfinite(ci[0]) else None,ci_high=float(ci[1]) if np.isfinite(ci[1]) else None,
        p_value=float(2*norm.sf(abs(effect/se))) if complete and np.isfinite(se) and se>0 else None,
        bootstrap_requested=bootstrap,bootstrap_success=len(valid),supported_donors=int(supported),omitted_markers=sorted(omitted),excluded_types=sorted(excluded_types),
        family='supporting_pseudotime',production_typing_modified=False,module_cache=str(module_output),
        uncertainty='Donor resampling reconstructs each fold diffusion graph and smooth; typing is fold-reconstructed and conditional on reviewed training identity references.')
    write_json(output/'result.json',result);return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True,type=Path);p.add_argument('--output-root',required=True,type=Path)
    p.add_argument('--reference-root',required=True,type=Path);p.add_argument('--outcome',action='append',required=True)
    p.add_argument('--compartment',choices=['core','peri','both'],default='both');p.add_argument('--workers',type=int,default=1);p.add_argument('--fit-workers',type=int,default=1);p.add_argument('--bootstrap',type=int,default=200);a=p.parse_args()
    out=offline_report_root(a.input_root,a.output_root);configure_runtime(out/'cache',a.workers);results=[]
    for outcome in a.outcome:
        directory=out/safe_name(outcome);directory.mkdir(parents=True,exist_ok=True)
        try:results.append(run(a.input_root,directory,a.reference_root,outcome,a.workers,a.fit_workers,a.bootstrap,a.compartment))
        except Exception as error:
            record=dict(outcome=outcome,status='unavailable',reason=str(error));write_json(directory/'result.json',record);results.append(record)
    table=pd.DataFrame(results)
    if 'p_value' in table:
        from statsmodels.stats.multitest import multipletests
        table['q_value']=np.nan;valid=table.p_value.notna()
        if valid.any():table.loc[valid,'q_value']=multipletests(table.loc[valid,'p_value'],method='fdr_bh')[1]
    table.to_csv(out/'supporting_pseudotime_results.csv',index=False)

if __name__=='__main__':main()
