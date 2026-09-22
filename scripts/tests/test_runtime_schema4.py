"""Canonical runtime summaries must agree with typed counts and mask areas."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'research'))
from build_runtime_schema4 import build
from common import MARKERS
from validate_schema4 import validate_count_arithmetic, validate_structural_areas
from finalize_candidate import startup_verify_artifact


def candidate(tmp_path):
    for folder in ['trajectory','typing','analysis']:
        (tmp_path/folder).mkdir()
    (tmp_path/'candidate.json').write_text(json.dumps({'schema_version':'4'}))
    (tmp_path/'typing/catalog.json').write_text(json.dumps({'base_types':[{'name':'Beta cell'},{'name':'CHGA-supported hormone-unresolved endocrine'}], 'states':[]}))
    pd.DataFrame({'case_id':['1','2']}).to_csv(tmp_path/'cohort.csv', index=False)
    pd.DataFrame({'case_id':['1','2'], 'donor_status':['ND','T1D'],'age':[20,30],'sex':['F','M']}).to_csv(tmp_path/'donor_metadata.csv', index=False)
    islets = pd.DataFrame({'case_id':['1','2'], 'islet_id':['1_a','2_b'], 'islet_key':['a','b'],'area_um2':[100.,200.],'peri_area_um2':[50.,80.]})
    islets.to_parquet(tmp_path/'analysis/islets.parquet', index=False)
    counts = []
    for row in islets.itertuples():
        for compartment in ['core','peri']:
            for phenotype,count in [('Beta cell',3), ('CHGA-supported hormone-unresolved endocrine',1)]:
                counts.append({'case_id':row.case_id,'islet_id':row.islet_id,'compartment':compartment,'phenotype':phenotype,'count':count,'validated':True,'n_measured':10})
    pd.DataFrame(counts).to_parquet(tmp_path/'trajectory/counts.parquet', index=False)
    pd.DataFrame({'case_id':['1','2'],'islet_id':['1_a','2_b'],'endocrine_count':[4,4],'validated':True}).to_parquet(tmp_path/'trajectory/endocrine_counts.parquet', index=False)
    marker=[]
    for row in islets.itertuples():
        for compartment in ['core','peri']:
            for protein in MARKERS:
                marker.append({'case_id':row.case_id,'islet_id':row.islet_id,'compartment':compartment,'population':'All cells','marker':protein,'n_cells':10,'positive_count':2.,'positive_available':protein!='CD45','mean_log1p':1.5})
    pd.DataFrame(marker).to_parquet(tmp_path/'trajectory/marker_summaries.parquet', index=False)
    structures=[]
    for row in islets.itertuples():
        for compartment,area in [('core',row.area_um2),('peri',row.peri_area_um2)]:
            for structure in ['Nerve','Capillary','Lymphatic']:
                structures.append({'case_id':row.case_id,'islet_id':row.islet_id,'compartment':compartment,'structure':structure,
                    'positive_area_um2':5.,'valid_area_um2':area,'available':structure!='Lymphatic'})
    pd.DataFrame(structures).to_parquet(tmp_path/'analysis/spatial_measurements.parquet', index=False)
    return tmp_path


def test_canonical_summary_uses_revised_counts_and_preserves_unavailable_masks(tmp_path):
    root=candidate(tmp_path)
    rows=build(root)
    assert rows['islets']==2
    summary=pd.read_parquet(root/'summary/islets.parquet')
    assert summary['prop_Beta cell'].tolist()==[.3,.3]
    assert summary['count_Beta cell'].tolist()==[3,3]
    assert summary['endocrine_count'].tolist()==[4,4]
    assert summary['gender'].tolist()==['F','M']
    assert summary['Multi_Pos'].isna().all()  # cannot recover overlap from marginal positives
    targets=pd.read_parquet(root/'summary/targets.parquet')
    assert targets.loc[targets['class']=='Lymphatic','area_um2'].isna().all()
    union=targets[(targets['class']=='Nerve') & (targets['type']=='islet_union')]
    np.testing.assert_allclose(union.area_density,[10/150,10/280])
    assert targets['count'].isna().all()
    markers=pd.read_parquet(root/'summary/markers.parquet')
    assert set(markers.marker)==set(MARKERS)
    assert markers.loc[markers.marker=='CD45','pos_count'].isna().all()
    assert set(pd.read_parquet(root/'summary/combined_pairs.parquet').stem)=={'Beta cell','CHGA-supported hormone-unresolved endocrine'}


def test_old_geometry_and_unknown_identity_fail_before_materialization(tmp_path):
    root=candidate(tmp_path)
    counts=pd.read_parquet(root/'trajectory/counts.parquet')
    counts.loc[0,'islet_id']='old_revision_islet'
    counts.to_parquet(root/'trajectory/counts.parquet', index=False)
    with pytest.raises(ValueError, match='revisions do not align'):
        build(root)
    assert not (root/'summary').exists()


def test_unvalidated_counts_never_become_zero(tmp_path):
    root=candidate(tmp_path)
    counts=pd.read_parquet(root/'trajectory/counts.parquet')
    counts.loc[counts.phenotype=='Beta cell','validated']=False
    counts.to_parquet(root/'trajectory/counts.parquet', index=False)
    build(root)
    summary=pd.read_parquet(root/'summary/islets.parquet')
    assert summary['count_Beta cell'].isna().all()
    assert summary['prop_Beta cell'].isna().all()


def test_independent_validation_checks_union_ratios_and_missing_structural_masks():
    counts=pd.DataFrame({'case_id':['1','1'],'islet_id':['1_a','1_a'],'compartment':['core','core'],
        'phenotype':['Beta cell','CHGA-supported hormone-unresolved endocrine'],'count':[3.,1.],
        'validated':[True,True],'n_typed':[4,4],'raw_ratio':[.75,.25],'plotted_value':np.log1p([.75,.25])})
    endocrine=pd.DataFrame({'islet_id':['1_a'],'endocrine_count':[4.]})
    validate_count_arithmetic(counts,endocrine,counts.phenotype.tolist())
    endocrine.endocrine_count=3
    with pytest.raises(ValueError,match='ratio mismatch'):validate_count_arithmetic(counts,endocrine,counts.phenotype.tolist())
    structure=pd.DataFrame({'case_id':['1'],'islet_id':['1_a'],'compartment':['core'],'structure':['Nerve'],
        'positive_area_um2':[np.nan],'valid_area_um2':[50.],'available':[False]})
    validate_structural_areas(structure)
    structure.positive_area_um2=0.
    with pytest.raises(ValueError,match='Unavailable mask'):validate_structural_areas(structure)


def test_schema4_startup_verification_skips_large_immutable_trees():
    for name in ['images', 'geometry', 'cells', 'islet_boundaries',
                 'provenance', 'typing', 'trajectory', 'spatial_tissue']:
        assert not startup_verify_artifact(name, '4')
    for name in ['cohort', 'comp', 'trajectory_coordinates',
                 'trajectory_feature_manifest', 'classifier_eligibility']:
        assert startup_verify_artifact(name, '4')
    # Other schemas retain startup checking for their ordinary trajectory tree.
    assert startup_verify_artifact('trajectory', '3')
