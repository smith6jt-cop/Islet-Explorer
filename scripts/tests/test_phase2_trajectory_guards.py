"""Regressions from the preservation-first donor-expansion handoff."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'research'))
import trajectory
from common import sha256
from stage_schema4_revision import stage
from revise_trajectory_reference_rank import clone
from validate_schema4 import validate_reference_percentiles


def features():
    rng=np.random.default_rng(12)
    frame=pd.DataFrame(rng.uniform(size=(12,3)),columns=[f'core__count__{t}' for t in 'ABC'])
    frame['islet_id']=[f'islet-{i}' for i in range(12)]
    frame['case_id']=['1']*6+['2']*6
    frame['area_um2']=100.
    frame['eligible']=True
    frame['trajectory_fit_supported']=True
    return frame


@pytest.mark.parametrize('missing', ['training', 'features'])
def test_fit_cannot_silently_remove_an_unsupported_donor(missing):
    frame=features()
    if missing=='training': frame.loc[frame.case_id=='2','trajectory_fit_supported']=False
    else:
        frame=frame.drop(columns='trajectory_fit_supported')
        frame.loc[frame.case_id=='2','core__count__A']=np.nan
    with pytest.raises(ValueError,match='refusing to remove cohort donors.*2'):
        trajectory.fit(frame,k=2)


def test_checkpoint_reuse_depends_on_training_mask_and_reference(tmp_path,monkeypatch):
    frame=features()
    reference=frame[['islet_id','case_id']].assign(pseudotime_raw=np.linspace(0,1,len(frame)))
    calls=[]
    def evaluate(frame,reference,family,name,params,training,columns):
        calls.append((family,name))
        return dict(family=family,replicate=name,success=True,absolute_spearman=1.)
    monkeypatch.setattr(trajectory,'_validation_one',evaluate)
    def run():
        return trajectory.validate(frame,reference,bootstrap=1,workers=1,checkpoint_root=tmp_path)
    run(); size=len(calls)
    assert size==8  # Two donor omissions, five sensitivities, one bootstrap.
    run(); assert len(calls)==size
    frame.loc[0,'trajectory_fit_supported']=False
    run(); assert len(calls)==2*size
    reference.loc[0,'pseudotime_raw']=.5
    run(); assert len(calls)==3*size
    reference.loc[0,'islet_id']='changed-reference-id'
    run(); assert len(calls)==4*size


def ranks():
    ids=['one','two','three','four']
    values=[.4,.4,.8,np.nan]
    # Construct the expected ordering independently of the implementation.
    order=sorted(range(3),key=lambda i:(values[i],hashlib.sha256(ids[i].encode()).hexdigest()))
    rank=np.full(4,np.nan)
    for position,i in enumerate(order): rank[i]=(position+.5)/3
    return pd.DataFrame(dict(islet_id=ids,pseudotime_raw=values,pseudotime_percentile=rank))


@pytest.mark.parametrize('swap', [(0,1),(0,2)])
def test_valid_grid_with_wrong_islet_assignments_is_rejected(swap):
    frame=ranks(); validate_reference_percentiles(frame,{'reference_percentile_version':2})
    a,b=swap
    frame.loc[[a,b],'pseudotime_percentile']=frame.loc[[b,a],'pseudotime_percentile'].to_numpy()
    with pytest.raises(ValueError,match='raw order or SHA-256'):
        validate_reference_percentiles(frame,{'reference_percentile_version':2})


@pytest.mark.parametrize('bad_ids', [[None,'b'],['a','a'],['a','']])
def test_reference_rank_rejects_missing_or_duplicate_ids_even_off_axis(bad_ids):
    with pytest.raises(ValueError,match='identifiers'):
        trajectory.fixed_reference_rank([.5,np.nan],bad_ids)


@pytest.mark.parametrize('artifact', ['model.npz','validation.json','graph.npz'])
def test_existing_partial_fit_cannot_be_overwritten(tmp_path,artifact):
    source=tmp_path/'source';source.mkdir()
    output=tmp_path/'output';output.mkdir();(output/artifact).write_bytes(b'prior artifact')
    with pytest.raises(FileExistsError,match='Fresh output'):
        trajectory.main(['--input-root',str(source),'--output-root',str(output)])
    assert (output/artifact).read_bytes()==b'prior artifact'


def test_main_rejects_donor_omitted_before_fitting(tmp_path):
    source=tmp_path/'source';(source/'trajectory').mkdir(parents=True)
    frame=features();path=source/'trajectory/composition_features.parquet';frame.to_parquet(path,index=False)
    manifest=dict(schema_version=4,representation='typed_composition',features_sha256=sha256(path),
                  columns=[c for c in frame if c.startswith('core__')])
    (source/'trajectory/feature_manifest.json').write_text(json.dumps(manifest))
    pd.DataFrame({'case_id':['1','2','3']}).to_csv(source/'cohort.csv',index=False)
    output=tmp_path/'fit'
    with pytest.raises(ValueError,match='incomplete cohort.*3'):
        trajectory.main(['--input-root',str(source),'--output-root',str(output)])
    assert not output.exists()


def test_isolated_main_preserves_all_rows_and_review_does_not_revoke_eligibility(tmp_path):
    source=tmp_path/'source';(source/'trajectory').mkdir(parents=True)
    rng=np.random.default_rng(31)
    frame=pd.DataFrame(rng.uniform(size=(40,3)),columns=[f'core__count__{t}' for t in 'ABC'])
    columns=list(frame)
    frame['islet_id']=[f'islet-{i}' for i in range(40)]
    frame['case_id']=np.repeat(['1','2','3','4'],10)
    frame['area_um2']=100.
    frame['eligible']=True
    frame['trajectory_feature_supported']=True
    frame['geometry_valid']=True
    frame['measured_cells']=10
    frame['trajectory_fit_supported']=True
    frame['exclusion_reason']=''
    frame.loc[[0,10,20,30],'measured_cells']=1
    frame.loc[[0,10,20,30],'trajectory_fit_supported']=False
    path=source/'trajectory/composition_features.parquet';frame.to_parquet(path,index=False)
    (source/'trajectory/feature_manifest.json').write_text(json.dumps(dict(
        schema_version=4,representation='typed_composition',features_sha256=sha256(path),
        columns=columns,primary_columns=columns)))
    pd.DataFrame({'case_id':['1','2','3','4']}).to_csv(source/'cohort.csv',index=False)
    geometry=source/'geometry.json';geometry.write_text('{}')
    pd.DataFrame({'case_id':['1','2','3','4'],'geometry':[str(geometry)]*4}).to_csv(source/'source_manifest.csv',index=False)
    review=frame[['islet_id','case_id']].assign(decision='unreviewed',geometry_sha256=sha256(geometry),review_note='Pending')
    review_path=source/'review.csv';review.to_csv(review_path,index=False)
    before={p.name:sha256(p) for p in source.rglob('*') if p.is_file()}
    output=tmp_path/'isolated-fit'
    trajectory.main(['--input-root',str(source),'--output-root',str(output),
                     '--islet-review',str(review_path),'--bootstrap','0'])
    coordinates=pd.read_parquet(output/'coordinates.parquet')
    assert len(coordinates)==40 and set(coordinates.case_id)=={'1','2','3','4'}
    assert coordinates.eligible.all() and coordinates.trajectory_fit_supported.sum()==36
    specification=json.loads((output/'specification.json').read_text())
    assert specification['identity_review_sha256']==sha256(review_path)
    assert (output/'identity_review.csv').is_file()
    validate_reference_percentiles(coordinates,specification)
    assert before=={p.name:sha256(p) for p in source.rglob('*') if p.is_file()}


@pytest.mark.parametrize('runner', [stage,clone])
def test_revision_derived_files_are_independent_and_source_cannot_contain_output(tmp_path,runner):
    source=tmp_path/'source';source.mkdir()
    (source/'candidate.json').write_text(json.dumps(dict(schema_version='4')))
    (source/'manifest.json').write_text(json.dumps(dict(schema_version='4',expected_images=21)))
    (source/'typing').mkdir();(source/'typing/cells.parquet').write_bytes(b'original typed cells')
    (source/'raw').mkdir();(source/'raw/cells.parquet').write_bytes(b'immutable raw cells')
    with pytest.raises(ValueError,match='separate from source'): runner(source,source/'child')
    target=tmp_path/'target';runner(source,target)
    (target/'typing/cells.parquet').write_bytes(b'revised typed cells')
    assert (source/'typing/cells.parquet').read_bytes()==b'original typed cells'
    assert (source/'raw/cells.parquet').samefile(target/'raw/cells.parquet')
    with pytest.raises(FileExistsError): runner(source,target)
