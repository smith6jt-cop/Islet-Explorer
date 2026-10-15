#!/usr/bin/env python3
"""Audit historical-sample/current-QC effects on strict RESTORE control maxima.

Diagnostic only: never changes the selected current calibration or gate calls.
The historical compartment sample is allowed only here, with raw Cell-mean and
fresh corrected UUID matching checked explicitly. Actual production fitting
continues to use native exports bound to the selected candidate.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import sha256,write_json
from restore_positivity import adapter,rv,lineage,quality_settings
from cell_typing import cell_eligibility
from threadpoolctl import threadpool_limits


def compare(root,screen,history,out,donors,markers):
    adapter.patch_upstream();spec=json.loads((screen/'specification.json').read_text());quality=spec['quality'];rows=[]
    baseline=dict(lineage.FIRST_PASS_PAIRS+lineage.LEVEL2_PAIRS)
    archived=pd.concat([pd.read_csv(history/name/'pair_evaluations.csv',dtype={'donor':str}) for name in ['screen_l2_clipped','screen_l2_unclipped']],ignore_index=True).set_index(['donor','target','reference'])
    first_pass=pd.read_csv(history/'screen_first_pass_valid/pair_evaluations.csv',dtype={'donor':str}).set_index(['donor','target','reference'])
    for cid in donors:
        native_path=history/'sample'/f'donor_id={cid}'/'data_0.parquet';historical=pd.read_parquet(native_path)
        corrected_path=root/'corrected'/f'{cid}-cells.parquet';raw_path=root/'raw'/f'{cid}-cells.parquet'
        corrected=pd.read_parquet(corrected_path);raw=pd.read_parquet(raw_path)
        raw.object_id=raw.object_id.str.removeprefix(cid+':');corrected.object_id=corrected.object_id.str.removeprefix(cid+':')
        usable=cell_eligibility(corrected,**quality);corrected['current_eligible']=usable.typing_eligible
        frame=historical[['object_id']+[c for c in historical if c.startswith(('Membrane__','Cell__'))]].copy()
        needed=sorted(set(markers)|{baseline[m] for m in markers})
        frame=frame.merge(raw[['object_id']+needed].rename(columns={m:'raw__'+m for m in needed}),on='object_id',validate='one_to_one')
        frame=frame.merge(corrected[['object_id','current_eligible']+needed].rename(columns={m:'redsea__'+m for m in needed}),on='object_id',validate='one_to_one')
        if len(frame)!=len(historical):raise ValueError('Historical diagnostic sample no longer matches current UUIDs exactly')
        for marker in needed:
            if not np.allclose(frame['Cell__'+marker],frame['raw__'+marker],rtol=2e-5,atol=2e-4,equal_nan=True):raise ValueError('Historical native measurement mismatch')
        fresh=pd.read_parquet(screen/'donors'/cid/'sample.parquet')
        for marker in markers:
            reference=baseline[marker];config=rv.PairValidationConfig(input_qc_upper_quantile=1. if marker in lineage.CLIP_UNCLIPPED_TARGETS else .995)
            for label,source in [('historical_sample_before_current_QC',frame),('historical_sample_with_current_QC',frame.loc[frame.current_eligible]),('fresh_sample_with_current_QC',fresh)]:
                record=dict(case_id=cid,marker=marker,reference=reference,sensitivity=label,n_cells=len(source),
                    historical_sample_sha256=sha256(native_path),corrected_sha256=sha256(corrected_path),raw_sha256=sha256(raw_path))
                try:
                    with threadpool_limits(limits=1):evaluation=rv.evaluate_locked_pair(source,cid,marker,reference,config=config,pair_review=rv.PairReview.ACCEPTED)
                    record.update(state=evaluation['state'],divisor=evaluation['candidate_full_maximum'],
                        reference_n=evaluation['reference_n'],target_n=evaluation['target_n'])
                    values=source['redsea__'+marker].to_numpy(float).copy();values[~evaluation['full_reference_control']]=-np.inf
                    ix=int(np.argmax(values));record['maximum_control_object_id']=source.iloc[ix].object_id if np.isfinite(values[ix]) else None
                    controls=values[np.isfinite(values)]
                    if len(controls):
                        historical_value=float(archived.loc[(cid,marker,reference),'candidate_full_maximum'])
                        q99,q995=np.quantile(controls,[.99,.995])
                        record.update(control_q99=float(q99),control_q995=float(q995),archived_candidate_full_maximum=historical_value,
                            archived_equals_true_maximum=bool(np.isclose(historical_value,controls.max(),rtol=1e-10,atol=1e-10)),
                            archived_equals_q995=bool(np.isclose(historical_value,q995,rtol=1e-10,atol=1e-10)),
                            rejected_patch_triggered=bool(q99>0 and controls.max()>5*q99))
                        if (cid,marker,reference) in first_pass.index:
                            value=float(first_pass.loc[(cid,marker,reference),'candidate_full_maximum'])
                            record.update(first_pass_candidate_full_maximum=value,
                                first_pass_equals_true_maximum=bool(np.isclose(value,controls.max(),rtol=1e-10,atol=1e-10)),
                                first_pass_equals_q995=bool(np.isclose(value,q995,rtol=1e-10,atol=1e-10)))
                except (ValueError,AssertionError) as exc:record.update(state='INPUT_QC_DEGENERATE',reason=str(exc),divisor=None)
                rows.append(record)
    out.mkdir(parents=True,exist_ok=True);pd.DataFrame(rows).to_csv(out/'sampling_qc_sensitivity.csv',index=False)
    write_json(out/'sampling_qc_sensitivity.json',dict(status='diagnostic_only',runner_sha256=sha256(__file__),
        primary_specification_sha256=sha256(screen/'specification.json'),donors=donors,markers=markers,
        historical_sampling='DuckDB hash(native UUID) %20 ==0; stored original sample',
        fresh_sampling='SHA256(20260920:case_id:UUID) first16hex %20 ==0; current QC',
        caveat='Fresh native Membrane+Cell export preserves double precision; historical sample cast means to float32. No donor gates or primary computational parameters changed.'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True);p.add_argument('--screen-root',type=Path,required=True)
    p.add_argument('--history-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--donor',action='append',required=True);p.add_argument('--marker',action='append',required=True);a=p.parse_args()
    compare(a.input_root,a.screen_root,a.history_root,a.output_root,a.donor,a.marker)


if __name__=='__main__':main()
