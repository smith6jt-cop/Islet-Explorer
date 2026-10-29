"""Filtered REDSEA must honor immutable ancestry and refuse stale/partial publication."""
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform/typing'))
import refresh_filtered_cells as F
from test_refresh_ready_cells import completed, fake_output


@pytest.fixture
def filtered(completed, monkeypatch):
    run, tag = completed
    base = F.RF.ready_inputs(run, tag)
    base['generation'] = 'uniform.new_detection_filter/1'
    base['filter_manifest_sha256'] = 'filtered-current'
    monkeypatch.setattr(F.NF, 'accepted_inputs', lambda run, tag: dict(base))
    return run, tag, base


def test_missing_generation_never_runs_redsea(tmp_path, monkeypatch):
    monkeypatch.setattr(F.NF, 'accepted_inputs', lambda *a: None)
    def forbidden(*a, **kw):
        pytest.fail('Incomplete filter generation must not enter REDSEA')
    monkeypatch.setattr(F.RS, 'process', forbidden)
    assert F.ready_inputs(tmp_path, 'A') is None
    with pytest.raises(ValueError, match='incomplete'):
        F.refresh_one(tmp_path, 'A')


def test_success_pins_filter_ancestry_and_reuses_current_cache(filtered, monkeypatch):
    run, tag, base = filtered
    calls = []
    def process(run, tag, force=False, dest_dir=None):
        assert force and dest_dir != F.CT.typing_dir(run) / 'redsea'
        calls.append(tag)
        fake_output(run, tag, dest_dir)
    monkeypatch.setattr(F.RS, 'process', process)
    assert F.refresh_one(run, tag)['status'] == 'refreshed'
    record = json.loads(F.provenance_path(run, tag).read_text())
    assert record['inputs']['filter_manifest_sha256'] == 'filtered-current'
    assert 'filtered_refresh_code' in record['inputs']['files']
    assert F.refresh_one(run, tag)['status'] == 'current'
    assert calls == [tag]


def test_filter_mutation_during_measurement_withholds_publication(filtered, monkeypatch):
    run, tag, base = filtered
    def process(run, tag, force=False, dest_dir=None):
        fake_output(run, tag, dest_dir)
        base['filter_manifest_sha256'] = 'changed-after-measurement'
    monkeypatch.setattr(F.RS, 'process', process)
    with pytest.raises(ValueError, match='inputs changed during REDSEA'):
        F.refresh_one(run, tag)
    assert not F.provenance_path(run, tag).exists()
    assert not (F.CT.typing_dir(run) / 'redsea' / f'{tag}.parquet').exists()


def test_adapter_code_change_requires_restart(filtered, monkeypatch, tmp_path):
    run, tag, _ = filtered
    path = tmp_path / 'adapter.py'
    path.write_text('version one')
    monkeypatch.setattr(F, 'CODE_PATHS', {'filtered_refresh_code': path})
    monkeypatch.setattr(F, 'LOADED_CODE', {'filtered_refresh_code': F.sha256(path)})
    assert F.ready_inputs(run, tag) is not None
    path.write_text('version two')
    with pytest.raises(ValueError, match='changed since import'):
        F.ready_inputs(run, tag)


def test_cohort_nulls_not_published_for_incomplete_filter(tmp_path, monkeypatch):
    monkeypatch.setattr(F, 'ready_inputs', lambda run, tag: None)
    assert F.publish_empty_nulls(tmp_path, ['A', 'B']) == ['A', 'B']
    assert not (tmp_path / 'typing/empty_area_nulls.csv').exists()
