#!/usr/bin/env python3
"""Fresh, explicit-cohort RESTORE positivity; distinct from hierarchical typing.

Uses the vendored locked Membrane+Cell NNMF scorer verbatim, with genuine native
QuPath means and newly corrected whole-cell means. Historical pair acceptance
is retained as biological reference evidence, never copied as a new threshold
or as acceptance of the six new donors. Calls require a current evidence-bound
review. All 31 output marker fields are nullable; four markers have no declared
exclusive negative-reference pair and remain unavailable through this method.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ.setdefault(name,'1')
REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'scripts/codex_phenocycler'))
# This adapter imports phenocycler before NumPy, preserving its BLAS policy.
import codex_restore_validation as adapter
import codex_lineage as lineage
from codex_panel import DONORS as ORIGINAL_DONORS
from phenocycler import restore_validation as rv
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from common import MARKERS,sha256,write_json
from cell_typing import cell_eligibility

SEED=20260920
NO_PAIR=['BActin','NaKATPase','BCatenin','ColIV']


def pair_name(target,reference):
    return target.replace('.','_')+'_from_'+reference.replace('.','_')


def historical_acceptance(history_root):
    """Reconcile actual first-pass CSV acceptance and the later user handoff."""
    first=history_root/'reviews/pair_reviews_first_pass.csv'
    handoff=REPO/'scripts/codex_phenocycler/docs/HANDOFF.md'
    table=pd.read_csv(first);accepted={(r.target,r.reference):dict(kind='accepted_first_pass_image_review',
        path=str(first.resolve()),sha256=sha256(first),evidence=r.evidence)
        for r in table.itertuples() if r.review=='ACCEPTED'}
    content=handoff.read_text()
    if 'Everything else: **ACCEPTED**.' not in content or 'PDPN 0.995 panel' not in content:
        raise ValueError('Recorded level-2 acceptance changed; reconcile the actual handoff before screening')
    # Round 1 accepted the declared baseline pair web, with explicit unresolved
    # PDPN re-review and sparse/channel cautions. The panCK hormone comparators
    # were dossier candidates, not replacement approved primary references.
    for pair in map(tuple,lineage.LEVEL2_PAIRS):
        if pair[0]=='PDPN':continue
        accepted.setdefault(pair,dict(kind='user_level2_round1_handoff',path=str(handoff),
            sha256=sha256(handoff),evidence='2026-07-24 image-review round 1: everything else accepted; named exceptions retained as current-review cautions.'))
    return accepted


def cautions(cid,target):
    reasons=[]
    if cid not in ORIGINAL_DONORS:reasons.append('new_donor_no_historical_image_acceptance')
    if target in lineage.PREVALENCE_COLLAPSE_MARKERS:reasons.append('sparse_target_mode_requires_image_supported_positive_population')
    if target in ['HLADR','INS']:reasons.append('historical_sparse_target_or_background_split_failure')
    if target=='PDPN':reasons.append('historical_clipped_PDPN_recheck_pending')
    if target=='CD20':reasons.append('historical_positive_calls_overcall_outside_CD45_supported_B_context')
    if cid=='6549' and target in ['KRT14','CD4','CD44','HLADR','PDPN']:
        reasons.append('historical_donor_channel_problem_requires_current_review')
    return reasons


def quality_settings(path):
    document=json.loads(path.read_text());config=document.get('config',document)
    return dict(minimum_cell_area_um2=float(config['minimum_cell_area_um2']),
                technical_exclusion_regions=config.get('technical_exclusion_regions',[]))


def sampled_cells(raw,corrected,cid,quality,modulus):
    """Select eligible raw UUIDs without marker prevalence or disease tuning."""
    raw=raw.copy();corrected=corrected.copy()
    for frame in [raw,corrected]:
        frame.case_id=frame.case_id.astype(str);frame.object_id=frame.object_id.astype(str)
        if frame.duplicated(['case_id','object_id']).any() or set(frame.case_id)!={cid}:
            raise ValueError('Duplicate UUID or wrong-cohort measurement input')
    if set(raw.object_id)!=set(corrected.object_id):raise ValueError('Fresh raw/corrected UUID universes differ')
    eligibility=cell_eligibility(corrected,**quality).typing_eligible
    good=corrected.loc[eligibility].copy()
    choose=good.object_id.map(lambda value:int(hashlib.sha256(f'{SEED}:{value}'.encode()).hexdigest()[:16],16)%modulus==0)
    good=good.loc[choose].sort_values('object_id').reset_index(drop=True)
    if len(good)<20:raise ValueError('Too few eligible deterministic RESTORE sample cells')
    return good,raw.set_index('object_id').loc[good.object_id].reset_index()


def join_compartments(native,raw,corrected,cid):
    native=native.copy();native.object_id=cid+':'+native.object_id.astype(str)
    if native.object_id.duplicated().any() or set(native.object_id)!=set(corrected.object_id):
        raise ValueError('Native compartment export must match the requested UUID sample exactly')
    frame=corrected[['case_id','object_id','x_source_px','y_source_px','Cell Area']+MARKERS].rename(columns={m:'redsea__'+m for m in MARKERS})
    frame=frame.merge(raw[['object_id']+MARKERS].rename(columns={m:'raw__'+m for m in MARKERS}),on='object_id',validate='one_to_one')
    frame=frame.merge(native,on='object_id',validate='one_to_one')
    for marker in MARKERS:
        column='Cell__'+marker
        if column not in frame or 'Membrane__'+marker not in frame:raise ValueError('Genuine native Membrane+Cell measurements are required')
        measured=frame['raw__'+marker].to_numpy(float);exported=frame[column].to_numpy(float)
        if not np.allclose(measured,exported,rtol=2e-5,atol=2e-4,equal_nan=True):
            raise ValueError(f'{cid}: native Cell mean does not match canonical raw {marker}; reconcile measurement/channel identity')
    return frame


def donor_screen(root,out,cid,quality_path,history_root,modulus,qupath):
    adapter.patch_upstream();quality=quality_settings(quality_path)
    folder=out/'donors'/cid;folder.mkdir(parents=True,exist_ok=True)
    raw_path=root/'raw'/f'{cid}-cells.parquet';corrected_path=root/'corrected'/f'{cid}-cells.parquet'
    native_manifest=root/'native_cells'/cid/'manifest.json';source=json.loads(native_manifest.read_text())
    source_path=Path(source['source_native'])
    if sha256(source_path)!=source['source_sha256']:raise ValueError('Native source hierarchy changed since verified export')
    accepted=historical_acceptance(history_root)
    fingerprint=dict(case_id=cid,raw_sha256=sha256(raw_path),corrected_sha256=sha256(corrected_path),
        native_sha256=source['source_sha256'],quality_sha256=sha256(quality_path),sample_modulus=modulus,
        adapter_sha256=sha256(__file__),upstream_sha256=sha256(Path(rv.__file__)),
        exporter_sha256=sha256(REPO/'scripts/cohort_expansion/export_restore_compartments.groovy'),
        history={str(pair):record for pair,record in accepted.items()})
    complete=folder/'complete.json'
    if complete.exists():
        if json.loads(complete.read_text())['fingerprint']!=fingerprint:raise ValueError('Screen inputs changed; use a fresh output root')
        return pd.read_csv(folder/'pair_evaluations.csv',dtype={'donor':str})
    corrected=pd.read_parquet(corrected_path);raw=pd.read_parquet(raw_path)
    chosen,raw_chosen=sampled_cells(raw,corrected,cid,quality,modulus)
    uuids=folder/'requested_uuids.txt';uuids.write_text(''.join(value.removeprefix(cid+':')+'\n' for value in chosen.object_id))
    markers=folder/'markers.json';write_json(markers,MARKERS)
    native_dir=folder/'native'
    if not (native_dir/'manifest.json').exists():
        command=[qupath,'script','--args',str(source_path),'--args',str(uuids),'--args',str(markers),
            '--args',str(native_dir),str(REPO/'scripts/cohort_expansion/export_restore_compartments.groovy')]
        env=os.environ.copy();env['JAVA_TOOL_OPTIONS']='-Xmx5G -XX:ActiveProcessorCount=1 -Djava.awt.headless=true'
        with (folder/'qupath_export.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
    proof=json.loads((native_dir/'manifest.json').read_text())
    if proof['requested_uuid_sha256']!=sha256(uuids) or proof['source_sha256']!=source['source_sha256'] or proof['table_sha256']!=sha256(native_dir/'compartments.tsv'):
        raise ValueError('Native compartment export hashes no longer match this run')
    frame=join_compartments(pd.read_csv(native_dir/'compartments.tsv',sep='\t',dtype={'object_id':str}),raw_chosen,chosen,cid)
    frame.to_parquet(folder/'sample.parquet',index=False)
    rows=[];thresholds=[];stability=[]
    for target,reference in lineage.ALL_SCREEN_PAIRS:
        pair=(target,reference);name=pair_name(*pair);checkpoint=folder/(name+'.json')
        if checkpoint.exists():
            saved=json.loads(checkpoint.read_text())
            if saved['fingerprint']!=fingerprint:raise ValueError('Pair checkpoint inputs changed')
            rows.append(saved['metrics']);thresholds.extend(saved['threshold_diagnostics']);stability.extend(saved['seed_diagnostics']);continue
        config=rv.PairValidationConfig(input_qc_upper_quantile=1. if target in lineage.CLIP_UNCLIPPED_TARGETS else .995)
        review=rv.PairReview.ACCEPTED if pair in accepted else rv.PairReview.UNREVIEWED
        pair_thresholds=[];pair_stability=[]
        try:
            with threadpool_limits(limits=1):evaluation=rv.evaluate_locked_pair(frame,cid,target,reference,config=config,pair_review=review)
            record=rv.locked_pair_metrics(evaluation)
            record={key:(None if isinstance(value,(float,np.floating)) and not np.isfinite(value) else value) for key,value in record.items()}
            for axis,key in [('target','target_threshold_grid'),('reference','reference_threshold_grid')]:
                for value in evaluation[key].to_dict('records'):pair_thresholds.append(dict(donor=cid,target=target,reference=reference,axis=axis,**value))
            for value in evaluation['stability_rows']:pair_stability.append(dict(donor=cid,target=target,reference=reference,**value))
            labels=frame[['object_id','x_source_px','y_source_px']].copy()
            for column in ['full_reference_control','full_target_population','full_threshold_positive','full_qc_retained']:
                labels[column]=evaluation[column]
            labels.to_parquet(folder/(name+'_sample_labels.parquet'),index=False)
        except (ValueError,AssertionError) as exc:
            record=dict(method_version=rv.METHOD_VERSION,donor=cid,target=target,reference=reference,
                state='INPUT_QC_DEGENERATE',state_reason=str(exc),candidate_full_maximum=None,canonical_divisor=None)
        record.update(current_review_status='awaiting_current_calibration_review',
            historical_pair_accepted=pair in accepted,historical_evidence=accepted.get(pair),
            cautions=cautions(cid,target),sample_cells=len(frame),input_qc_upper_quantile=config.input_qc_upper_quantile)
        for records in [pair_thresholds,pair_stability]:
            for item in records:
                for key,value in item.items():
                    if isinstance(value,(float,np.floating)) and not np.isfinite(value):item[key]=None
        write_json(checkpoint,dict(fingerprint=fingerprint,metrics=record,threshold_diagnostics=pair_thresholds,seed_diagnostics=pair_stability))
        thresholds.extend(pair_thresholds);stability.extend(pair_stability)
        rows.append(record)
    result=pd.DataFrame(rows);result.to_csv(folder/'pair_evaluations.csv',index=False)
    pd.DataFrame(thresholds).to_csv(folder/'threshold_sampling.csv',index=False)
    pd.DataFrame(stability).to_csv(folder/'seed_stability.csv',index=False)
    write_json(complete,dict(fingerprint=fingerprint,sample_sha256=sha256(folder/'sample.parquet'),
        current_calibration_review='pending',sample_cells=len(frame),locked_method_version=rv.METHOD_VERSION,
        divisor_policy='Full maximum corrected target intensity among inferred mutually exclusive reference controls. No quantile substitution or automatic donor borrowing.'))
    print(cid,result.state.value_counts().to_dict(),flush=True)
    return result


def discover(root,out,quality_path,history_root,modulus,qupath,workers,selected=None):
    cohort=pd.read_csv(root/'cohort.csv',dtype={'case_id':str});donors=sorted(cohort.case_id)
    if cohort.case_id.duplicated().any():raise ValueError('Duplicate cohort donors')
    if selected:
        if set(selected)-set(donors):raise ValueError('Requested donor outside explicit cohort')
        donors=[cid for cid in donors if cid in selected]
    out.mkdir(parents=True,exist_ok=True);frames=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs=[pool.submit(donor_screen,root,out,cid,quality_path,history_root,modulus,qupath) for cid in donors]
        for job in as_completed(jobs):frames.append(job.result())
    result=pd.concat(frames,ignore_index=True);result.to_csv(out/'pair_evaluations.csv',index=False)
    write_json(out/'specification.json',dict(schema_version=1,method='fresh_native_Membrane_Cell_locked_RESTORE_v8',
        method_version=rv.METHOD_VERSION,cohort_sha256=sha256(root/'cohort.csv'),donors=donors,
        cohort_donors=sorted(cohort.case_id),complete_cohort=set(donors)==set(cohort.case_id),seed=SEED,sample_modulus=modulus,
        quality_specification=str(quality_path.resolve()),quality_sha256=sha256(quality_path),quality=quality_settings(quality_path),
        runner_sha256=sha256(__file__),no_declared_negative_reference=NO_PAIR,
        historical_acceptance_scope='Original15 biological pair choices, reconciled with user level2 handoff; no old numeric threshold copied.',
        positivity_call='fresh corrected whole-cell value > reviewed donor-specific canonical divisor',
        state='fresh_computation_complete_current_calibration_review_pending',typing_role='none; independent marker summaries and count endpoints only'))
    template=dict(specification_sha256=sha256(out/'specification.json'),reviewer=None,donors={})
    for cid in donors:
        template['donors'][cid]={m:dict(validated=False,reference=None,evaluation_sha256=None,evidence=[],reason='Current calibration review pending' if m not in NO_PAIR else 'No declared exclusive reference population') for m in MARKERS}
    write_json(out/'review_template.json',template)


def apply(root,out,review_path,destination):
    spec=json.loads((out/'specification.json').read_text());review=json.loads(review_path.read_text())
    if sha256(out/'specification.json')!=review['specification_sha256'] or not review.get('reviewer'):
        raise ValueError('Current reviewer and exact specification binding required')
    donors=pd.read_csv(root/'cohort.csv',dtype=str).case_id.tolist()
    if not spec['complete_cohort'] or set(review['donors'])!=set(donors):raise ValueError('Production positivity must explicitly account for the complete manifest cohort')
    if sha256(Path(spec['quality_specification']))!=spec['quality_sha256']:raise ValueError('QC specification changed')
    if sha256(root/'cohort.csv')!=spec['cohort_sha256']:raise ValueError('Cohort changed after calibration')
    destination.mkdir(parents=True,exist_ok=True);rows=[];outputs={};sources={}
    for cid in donors:
        decisions=review['donors'][cid]
        if set(decisions)!=set(MARKERS):raise ValueError('Every marker needs an explicit review disposition')
        checkpoint=json.loads((out/'donors'/cid/'complete.json').read_text())
        path=root/'corrected'/f'{cid}-cells.parquet'
        if sha256(path)!=checkpoint['fingerprint']['corrected_sha256']:raise ValueError('Corrected cells changed since calibration')
        sources[cid]=dict(path=str(path.resolve()),sha256=checkpoint['fingerprint']['corrected_sha256'])
        cells=pd.read_parquet(path);eligible=cell_eligibility(cells,**spec['quality']).typing_eligible
        calls=cells[['case_id','object_id']].copy()
        for marker,decision in decisions.items():
            values=pd.Series(pd.NA,index=cells.index,dtype='boolean')
            if not decision.get('reason'):raise ValueError('Review explanation required')
            threshold=None
            if decision['validated']:
                if marker in NO_PAIR:raise ValueError('Cannot validate an undeclared negative-reference pair with RESTORE')
                if not decision['evidence']:raise ValueError('Current calibration needs actual image evidence')
                for evidence in decision['evidence']:
                    if sha256(evidence['path'])!=evidence['sha256']:raise ValueError('Current image evidence changed')
                pair_path=out/'donors'/cid/(pair_name(marker,decision['reference'])+'.json')
                if sha256(pair_path)!=decision['evaluation_sha256']:raise ValueError('Reviewed divisor belongs to another computation')
                evaluation=json.loads(pair_path.read_text())['metrics']
                if evaluation['state'] not in ['VALID','PAIR_REVIEW_PENDING']:raise ValueError('Review cannot overrule failed computational prerequisites')
                threshold=evaluation.get('candidate_full_maximum')
                if threshold is None or not np.isfinite(threshold) or threshold<=0:raise ValueError('Missing positive fresh control maximum')
                values.loc[eligible]=cells.loc[eligible,marker].gt(threshold)
            calls['positive__'+marker]=values
            rows.append(dict(case_id=cid,marker=marker,validated=decision['validated'],threshold=threshold,
                             reference=decision.get('reference'),reason=decision['reason'],eligible_cells=int(eligible.sum())))
        target=destination/f'{cid}.parquet'
        if target.exists():raise FileExistsError('Preserve previous positivity calls; choose a fresh destination')
        calls.to_parquet(target,index=False)
        outputs[cid]=dict(path=str(target.resolve()),sha256=sha256(target),rows=len(calls),eligible_cells=int(eligible.sum()))
    pd.DataFrame(rows).to_csv(destination/'marker_availability.csv',index=False)
    write_json(destination/'manifest.json',dict(status='current_reviewed_RESTORE_calls',schema_version=1,
        review_path=str(review_path.resolve()),review_sha256=sha256(review_path),source_specification_path=str((out/'specification.json').resolve()),source_specification_sha256=sha256(out/'specification.json'),
        cohort_manifest_path=str((root/'cohort.csv').resolve()),cohort_manifest_sha256=spec['cohort_sha256'],
        quality_specification_path=spec['quality_specification'],quality_sha256=spec['quality_sha256'],quality=spec['quality'],
        corrected_sources=sources,outputs=outputs,calls_sha256={cid:value["sha256"] for cid,value in outputs.items()},corrected_inputs={cid:value["sha256"] for cid,value in sources.items()},apply_runner_sha256=sha256(__file__),
        donors=donors,no_declared_negative_reference=NO_PAIR,typing_role='none',disease_used=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['discover','apply'])
    parser.add_argument('--input-root',type=Path,required=True);parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--typing-spec',type=Path);parser.add_argument('--history-root',type=Path)
    parser.add_argument('--workers',type=int,default=2);parser.add_argument('--sample-modulus',type=int,default=20)
    parser.add_argument('--donor',action='append');parser.add_argument('--qupath',default='/home/smith6jt/QuPath/bin/QuPath')
    parser.add_argument('--review',type=Path);parser.add_argument('--calls-root',type=Path)
    args=parser.parse_args();root=args.input_root.resolve(strict=True);out=args.output_root.resolve()
    if args.command=='discover':
        if not args.typing_spec or not args.history_root:parser.error('discover requires --typing-spec and --history-root')
        if args.sample_modulus<1:parser.error('sample-modulus must be positive')
        discover(root,out,args.typing_spec.resolve(strict=True),args.history_root.resolve(strict=True),args.sample_modulus,args.qupath,args.workers,args.donor)
    else:
        if not args.review or not args.calls_root:parser.error('apply requires --review and --calls-root')
        apply(root,out,args.review,args.calls_root)


if __name__=='__main__':main()
