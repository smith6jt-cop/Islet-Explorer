#!/usr/bin/env python3
"""Independent candidate checks, never reuse the frozen release's evidence."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import MARKERS, RETIRED, write_json, require_cohort, sha256

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--report-root',type=Path);a=p.parse_args();root=Path(a.input_root)
    cohort=pd.read_csv(root/'cohort.csv',dtype=str);report={}
    assert len(cohort)==21 and cohort.donor_status.value_counts().to_dict()=={'T1D':7,'ND':7,'Aab+':7}
    assert cohort.loc[cohort.donor_source=='HPAP','case_id'].tolist()==['78']
    images=sorted(p.name for p in (root/'images').glob('*.ome.tiff'))
    assert images==sorted(Path(p).name for p in cohort.image)
    cells=unmeasured=0
    for cid in cohort.case_id:
        c=pd.read_parquet(root/'raw'/f'{cid}-cells.parquet',columns=['case_id','object_id','marker_measurements_complete']+MARKERS)
        assert c.object_id.is_unique and set(c.case_id)=={cid}
        finite=np.isfinite(c[MARKERS]);measured=finite.all(axis=1)
        assert (measured==c.marker_measurements_complete).all()
        assert ((finite.sum(axis=1)==0)|measured).all() and (c.loc[measured,MARKERS]>=0).all().all()
        cells+=len(c);unmeasured+=int((~measured).sum())
    assert cells==3592722 and unmeasured==201
    report.update(donors=21,raw_objects=cells,wholly_unmeasured=unmeasured,shared_markers=31)
    r=pd.read_parquet(root/'analysis/regions.parquet');require_cohort(r,cohort.case_id)
    assert not r.duplicated(['islet_id','compartment','band_um']).any()
    assert set(r.band_um)=={10,20,30}
    for s in ['nerve','capillary','lymphatic']:
        good=r.valid_area_um2>0;actual=r.loc[good,s+'_area_um2']/r.loc[good,'valid_area_um2']
        assert np.allclose(actual,r.loc[good,s+'_area_fraction'])
        assert actual.min()>=-1e-8 and actual.max()<=1+1e-8
    f=pd.read_csv(root/'trajectory/features.csv',dtype={'case_id':str});c=pd.read_csv(root/'trajectory/coordinates.csv',dtype={'case_id':str})
    assert f.islet_id.is_unique and set(c.islet_id)==set(f.islet_id)
    lookup=json.loads((root/'cell_lookup.json').read_text())
    assert {x['case_id'] for x in lookup['donors']}==set(cohort.case_id)
    cell_files={p.stem:p for p in (root/'cells').glob('*.csv')}
    assert set(f.loc[f.eligible,'islet_id']).issubset(cell_files)
    assert {s.split('_')[0] for s in cell_files}==set(cohort.case_id)
    for cid in cohort.case_id:
        path=next(p for s,p in cell_files.items() if s.startswith(cid+'_'))
        columns=pd.read_csv(path,nrows=0).columns
        assert set(MARKERS+['object_id','X_centroid','Y_centroid','marker_measurements_complete','cell_region']).issubset(columns)
    report.update(cell_drilldown_files=len(cell_files),cell_drilldown_donors=len(lookup['donors']))
    assert (f.loc[f.eligible,'measured_cells']>=10).all()
    features=[x for x in f if x.startswith(('core__','peri__'))]
    assert len(features)==62 and np.isfinite(f.loc[f.eligible,features]).all().all()
    primary=pd.read_json(root/'trajectory/specification.json',typ='series')
    assert not {'INS','GCG','SST','CHGA'} & {s.split('__')[1] for s in primary['columns']}
    counts=pd.read_csv(root/'trajectory/counts.csv');assert not set(counts.phenotype)&set(RETIRED)
    assert not counts.loc[~counts.validated,'count'].notna().any()
    validation=json.loads((root/'trajectory/validation.json').read_text())
    assert validation['bootstrap_completed']==200 and validation['lodo_completed']==21
    report.update(islets=len(f),eligible=int(f.eligible.sum()),geometry_area_fraction_checked=len(r),
                  trajectory_status=validation['status'],phenotype_uncertainty_preserved=True)
    if (root/'trajectory/endocrine_counts.csv').exists():
        denom=pd.read_csv(root/'trajectory/endocrine_counts.csv',dtype={'case_id':str})
        require_cohort(denom,cohort.case_id)
        assert denom.islet_id.is_unique and set(denom.islet_id)==set(f.islet_id)
        assert denom.loc[~denom.validated,'endocrine_count'].isna().all()
        assert (denom.loc[denom.validated,'endocrine_count']>=0).all()
        assert not counts.duplicated(['islet_id','phenotype','compartment']).any()
        assert set(counts.islet_id)==set(f.islet_id)
        panel=json.loads((root/'trajectory/marker_panel.json').read_text())
        assert panel['markers']==MARKERS and set('core__'+m for m in MARKERS).issubset(f)
        assert not (root/'research').exists()
        report.update(research_not_in_runtime=True,validated_endocrine_donors=int(denom.loc[denom.validated,'case_id'].nunique()))
    if (root/'revision.json').exists():
        revision=json.loads((root/'revision.json').read_text());previous=Path(revision['source_root'])
        files=['features.csv','coordinates.csv','edges.csv','graph.npz','landmarks.csv','model.npz','specification.json','validation.json','validation_replicates.csv']
        for file in files:assert sha256(root/'trajectory'/file)==sha256(previous/'trajectory'/file),file
        report['fitted_trajectory_artifacts_unchanged']=files
    if (root/'trajectory/quality_review.json').exists():
        quality=json.loads((root/'trajectory/quality_review.json').read_text())
        report['scientific_status']=quality['status'];report['scientific_acceptance']=False
    if a.report_root:
        assert not a.report_root.resolve().is_relative_to(root.resolve())
        results=pd.read_csv(a.report_root/'results.csv')
        assert set(results.id)=={f'H{i}' for i in range(1,11)}
        assert results.loc[results.status=='unavailable','p_value'].isna().all()
        report.update(estimated_primary_tests=int(results.p_value.notna().sum()),unavailable_primary_tests=int(results.p_value.isna().sum()))
    report.update(passed=True,validation_script_sha256=sha256(__file__))
    write_json(root/'candidate_validation.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
