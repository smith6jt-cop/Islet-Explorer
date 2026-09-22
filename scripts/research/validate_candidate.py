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
    identity=json.loads((root/'candidate.json').read_text())
    if str(identity.get('schema_version'))=='4':
        from validate_schema4 import validate
        print(json.dumps(validate(root,a.report_root),indent=2))
        return
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
    for name in ['features','coordinates','counts','endocrine_counts','edges']:
        parquet=root/'trajectory'/f'{name}.parquet'
        if parquet.exists():
            csv=pd.read_csv(parquet.with_suffix('.csv'),dtype={'case_id':str})
            pd.testing.assert_frame_equal(csv,pd.read_parquet(parquet),check_dtype=False)
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
    if validation['status']=='awaiting_reviewed_refit':
        assert validation['coordinates_available'] is False
        assert validation['bootstrap_completed']==0 and validation['lodo_completed']==0
        assert c[['pseudotime_raw','pseudotime_percentile']+['DC'+str(i) for i in range(1,11)]].isna().all().all()
        assert not c.supported.any() and not c.landmark.any() and (c.component==-1).all()
        assert pd.read_csv(root/'trajectory/edges.csv').empty
        assert not any((root/'trajectory'/name).exists() for name in ['model.npz','graph.npz','landmarks.csv','validation_replicates.csv'])
        assert primary['features_sha256']==sha256(root/'trajectory/features.csv')
        assert primary['cohort_sha256']==sha256(root/'cohort.csv')
        assert primary['islet_review_sha256']==sha256(root/'trajectory/islet_review.csv')
        assert primary['previous_model_reused'] is False
        assert primary['source_geometry_sha256']=={cid:sha256(root/'gson'/f'{int(cid):04d}.geojson.gz') for cid in cohort.case_id}
    else:assert validation['bootstrap_completed']==200 and validation['lodo_completed']==21
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
        if revision.get('kind')=='approved_boundary_integration':
            assert sha256(previous/'manifest.json')==revision['source_manifest_sha256']
            assert validation['status']=='awaiting_reviewed_refit'
            hpap=f[f.case_id=='78'];assert len(hpap)==1431 and hpap.islet_key.str.startswith('AdaptiveIslet_78_').all()
            assert len(f)==8568 and not list((root/'cells').glob('78_Islet_*.csv'))
            assert {s for s in cell_files if s.startswith('78_')}==set(hpap.islet_id)
            native=pd.read_parquet(root/'analysis/78-memberships.parquet')
            expected=set(map(tuple,native[['islet_id','object_id']].values));actual=set()
            assert len(expected)==9569
            for row in hpap.itertuples():
                frame=pd.read_csv(cell_files[row.islet_id])
                core=frame.loc[frame.cell_region=='core']
                assert core.object_id.is_unique
                assert int(core.marker_measurements_complete.sum())==row.measured_cells
                actual.update((row.islet_id,oid) for oid in core.object_id)
            assert actual==expected
            import sqlite3
            with sqlite3.connect(root/'geometry/0078.gpkg') as con:
                indexed={r[0] for r in con.execute('SELECT name FROM segmentation WHERE name LIKE "AdaptiveIslet_78_%"')}
                assert set(hpap.islet_key).issubset(indexed)
                assert con.execute('SELECT count(*) FROM segmentation WHERE name LIKE "Islet_%"').fetchone()[0]==0
            for path in (root/'summary').glob('*.parquet'):
                frame=pd.read_parquet(path);old=pd.read_parquet(previous/'summary'/path.name)
                pd.testing.assert_frame_equal(frame.loc[frame['Case ID']!=78].reset_index(drop=True),old.loc[old['Case ID']!=78].reset_index(drop=True))
                if not frame.empty:assert set(frame.loc[frame['Case ID']==78,'islet_key'])==set(hpap.islet_key)
            old=pd.read_csv(previous/'trajectory/features.csv',dtype={'case_id':str})
            columns=old.columns
            pd.testing.assert_frame_equal(f.loc[f.case_id!='78',columns].reset_index(drop=True),old.loc[old.case_id!='78'].reset_index(drop=True))
            review=pd.read_csv(root/'trajectory/boundary_review.csv',dtype={'case_id':str})
            assert review.loc[review.approved,'case_id'].tolist()==['78'] and review.pending_new_candidates.sum()==2030
            assert not (root/'master_results.xlsx').exists()
            if (root/'islet_boundaries').exists():
                import gzip
                outlines={p.stem:p for p in (root/'islet_boundaries').glob('*.geojson')}
                assert set(outlines)==set(f.islet_id)
                with gzip.open(root/'gson/0078.geojson.gz','rt') as stream:original=json.load(stream)
                native_features={v['properties'].get('name'):v for v in original['features'] if v['properties'].get('classification',{}).get('name')=='Islet'}
                for row in hpap.itertuples():
                    downloaded=json.loads(outlines[row.islet_id].read_text())
                    assert downloaded['features']==[native_features[row.islet_key]]
                report.update(single_islet_outlines=len(outlines),exact_hpap_download_vertices_verified=True)
            report.update(corrected_hpap_islets=1431,native_core_memberships_verified=len(actual),
                          previous_trajectory_withdrawn=True,other_20_donor_measurements_preserved=True)
        else:
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
