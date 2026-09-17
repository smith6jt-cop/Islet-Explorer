#!/usr/bin/env python3
"""Endocrine marker profiles using the literal QuPath hormone thresholds.

This is an operational hormone-positive union, independent of the historical
phenotype parent gates. It is not substituted for validated subtype assignments
in the trajectory ratio denominator.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from common import MARKERS, write_json, offline_report_root

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--report-root',required=True);a=p.parse_args();root=Path(a.input_root);out=offline_report_root(root,a.report_root)
    donors=pd.read_csv(root/'donor_metadata.csv',dtype={'case_id':str});rows=[];audit=[]
    for cid in donors.case_id:
        cells=pd.read_parquet(root/'raw'/f'{cid}-cells.parquet').set_index('object_id')
        composition=pd.read_csv(root/'regional'/f'islet_composition_{int(cid):04d}.tsv',sep='\t')
        thresholds=composition[['Ins_thr','Glu_thr','Stt_thr']].drop_duplicates()
        if len(thresholds)!=1:raise ValueError('Expected one explicitly exported hormone gate per donor')
        union=np.any(cells[['INS','GCG','SST']].to_numpy()>thresholds.iloc[0].to_numpy(),axis=1)&cells.marker_measurements_complete.to_numpy()
        logs=np.log1p(cells[MARKERS].to_numpy())
        membership=pd.read_parquet(root/'analysis'/f'{cid}-memberships.parquet')
        expected=(composition.set_index('name').cells_total-composition.set_index('name').Triple_Neg)
        for islet,group in membership.groupby('islet_id'):
            ix=cells.index.get_indexer(group.object_id)
            if np.any(ix<0):raise ValueError('Cell membership points outside selected donor')
            valid=ix[union[ix]];name=islet[len(cid)+1:]
            source_count=int(expected.loc[name]);delta=len(valid)-source_count
            audit.append(dict(islet_id=islet,case_id=cid,raw_threshold_union=len(valid),exported_union=source_count,count_difference=delta))
            record=dict(islet_id=islet,case_id=cid,endocrine_measured_cells=len(valid),count_difference=delta,
                        profile_eligible=len(valid)>=5 and delta==0,
                        median_cell_area=float(pd.to_numeric(cells.iloc[valid]['Cell Area']).median()) if len(valid) else np.nan)
            for m,v in zip(MARKERS,np.median(logs[valid],axis=0) if len(valid) else np.full(31,np.nan)):record[m]=v
            rows.append(record)
        print(cid,'endocrine profiles',flush=True)
    pd.DataFrame(rows).to_parquet(root/'analysis/endocrine_profiles.parquet',index=False)
    pd.DataFrame(audit).to_csv(out/'endocrine_threshold_audit.csv',index=False)
    write_json(out/'endocrine_profile_specification.json',dict(selection='INS > Ins_thr OR GCG > Glu_thr OR SST > Stt_thr; measured cells only',
               minimum_endocrine_cells=5,primary_count_validation='exact agreement with independent QuPath exported hormone-positive union',
               markers=MARKERS,phenotype_subtype_denominator=False,
               limitation='Baseline centroid rounding and rare float32 threshold ties can cause mismatches; excluded from primary endocrine-profile analysis and retained in audit.'))

if __name__=='__main__':main()
