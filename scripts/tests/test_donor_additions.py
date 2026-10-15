"""Raw staging cannot relabel measurements, drop cells, or mutate sources."""
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'research'))
import stage_donor_additions as stage
from common import MARKERS, require_fresh_output, sha256


def raw(path, case='78'):
    obs = pd.DataFrame({'imageid': [case]*3, 'Classification': ['INS_pos', '', 'Multi_Pos']},
                       index=['uuid1', 'uuid2', 'uuid3'])
    values = np.ones((3, len(MARKERS)), dtype=np.float32)
    values[1] = np.nan
    a = ad.AnnData(values, obs=obs, var=pd.DataFrame(index=MARKERS))
    a.write_h5ad(path)
    return a


def test_raw_snapshot_preserves_missingness_and_archival_labels(tmp_path):
    source, target = tmp_path/'raw.h5ad', tmp_path/'copy.h5ad'
    raw(source)
    before = sha256(source)
    stage.copy_verified(source, target, before)
    census = stage.inspect_raw(target, '78')
    assert census['raw_cells'] == 3 and census['wholly_unmeasured_cells'] == 1
    a = ad.read_h5ad(target)
    assert 'phenotype' not in a.obs and not a.obsm
    assert a.obs.Classification.iloc[0] == 'INS_pos'
    assert a.obs_names.tolist() == ['uuid1', 'uuid2', 'uuid3']
    a.X[0] = 0
    a.write_h5ad(target)
    assert sha256(source) == before


def test_stale_import_hash_creates_no_copy(tmp_path):
    source, target = tmp_path/'raw.h5ad', tmp_path/'copy.h5ad'
    raw(source)
    with pytest.raises(ValueError, match='Source hash'):
        stage.copy_verified(source, target, 'stale')
    assert not target.exists()


@pytest.mark.parametrize('defect', ['donor', 'duplicate', 'phenotype', 'scvi', 'infinity'])
def test_reject_incompatible_raw_sources(tmp_path, defect):
    path = tmp_path/'source.h5ad'; a = raw(path)
    if defect == 'donor': a.obs['imageid'] = '6090'
    if defect == 'duplicate': a.obs_names = ['uuid1']*3
    if defect == 'phenotype': a.obs['phenotype'] = 'Beta cell'
    if defect == 'scvi': a.obsm['X_scVI'] = np.zeros((3, 2))
    if defect == 'infinity': a.X[0, 0] = np.inf
    a.write_h5ad(path)
    with pytest.raises(ValueError): stage.inspect_raw(path, '78')


def test_fresh_output_rejects_existing_recursive_and_broken_links(tmp_path):
    source = tmp_path/'source'; source.mkdir()
    broken = tmp_path/'broken'; broken.symlink_to(tmp_path/'absent')
    for target in [source, broken]:
        with pytest.raises(FileExistsError): require_fresh_output(target, [source])
    with pytest.raises(ValueError): require_fresh_output(source/'child', [source])
    alias = tmp_path/'alias'; alias.symlink_to(source, target_is_directory=True)
    with pytest.raises(ValueError): require_fresh_output(alias/'child', [source])
    assert not (source/'child').exists()


def test_full_staging_uses_only_six_additions_and_independent_originals(tmp_path, monkeypatch):
    baseline = tmp_path/'original'; baseline.mkdir()
    original = ad.AnnData(np.ones((15, 1)), obs=pd.DataFrame(
        {'imageid': [str(i) for i in range(1, 16)]}, index=[str(i) for i in range(15)]))
    original.write_h5ad(baseline/'adata_ins_root.h5ad')
    for name in stage.BASELINE_FILES[1:]: (baseline/name).write_text('frozen')
    rows = []
    for cid in sorted(stage.ADDITIONS):
        folder = tmp_path/'inputs'/cid; folder.mkdir(parents=True)
        path = folder/'cells.h5ad'; raw(path, cid)
        path.with_suffix('.json').write_text(json.dumps({'output_sha256': sha256(path)}))
        (folder/'geometry.gz').write_bytes(b'geometry')
        (folder/'metadata.csv').write_text('metadata')
        rows.append(dict(case_id=cid, cells=str(path), cell_slot='X', geometry=str(folder/'geometry.gz'),
                         metadata=str(folder/'metadata.csv'), donor_source='HPAP' if cid=='78' else 'nPOD',
                         donor_status='T1D'))
    rows += [dict(case_id=str(i)) for i in range(1, 16)]
    cohort = pd.DataFrame(rows); cohort_path=tmp_path/'cohort.csv'; cohort.to_csv(cohort_path, index=False)
    monkeypatch.setattr(stage, 'read_cohort', lambda *args: cohort)
    monkeypatch.setattr(stage, 'EXPECTED_TOTALS', dict(raw_cells=18, fully_measured_cells=12,
                        wholly_unmeasured_cells=6, partially_measured_cells=0))
    output = tmp_path/'bundle'
    result = stage.stage(cohort_path, tmp_path, baseline, output)
    assert result['totals']['raw_cells'] == 18 and result['app_ready'] is False
    assert {p.stem for p in (output/'additions/raw').iterdir()} == stage.ADDITIONS
    for name in stage.BASELINE_FILES:
        assert sha256(output/name) == sha256(baseline/name)
        assert not (output/name).samefile(baseline/name)
    assert not (output/'candidate.json').exists()
    with pytest.raises(FileExistsError): stage.stage(cohort_path, tmp_path, baseline, output)
