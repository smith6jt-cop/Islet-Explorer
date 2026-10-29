"""The overlap watcher must not measure unfinished augmentation or reuse copied REDSEA."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import sys
import threading

import pandas as pd
import pytest

TYPING = Path(__file__).resolve().parents[1] / 'uniform' / 'typing'
sys.path.insert(0, str(TYPING))
import refresh_ready_cells as RF


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def completed(tmp_path, monkeypatch):
    run = tmp_path / 'run'
    tag = 'A'
    image = tmp_path / 'image.tif'
    image.write_bytes(b'immutable original image')
    write_json(run / 'project_manifest.json', {'entries': [
        {'image_name': tag, 'entry_id': 1, 'image_path': str(image)}]})
    write_json(run / 'qupath/data/1/server.json', {'metadata': {'channels': [{'name': 'DAPI'}]}})
    cells = run / 'cells' / f'{tag}.parquet'
    cells.parent.mkdir(parents=True)
    pd.DataFrame({'Object ID': ['old', 'new'], 'Cell: DAPI: Mean': [5., 7.]}).to_parquet(cells, index=False)
    geo = run / 'cell_geojson' / f'{tag}.geojson'
    tissue = run / 'inputs/tissue' / f'{tag}.geojson'
    report = run / 'augmentation' / tag / 'verification.json'
    write_json(geo, {'features': [{'id': 'old'}, {'id': 'new'}]})
    write_json(tissue, {'features': []})
    write_json(report, {'original_cells': 1, 'new_cells': 1, 'total_cells': 2})
    state = {'schema': 'uniform.tissue_augmentation/1', 'inputs': {'outline_sha256': RF.sha256(tissue)},
             'last_qpdata_sha256': 'verified-qpdata',
             'steps': {s: {'fingerprint': 'completed', 'outputs': {}} for s in RF.REQUIRED_STEPS}}
    for step, path in [('augment_detect', tissue), ('augment_C_cells', cells),
                       ('augment_C_geojson', geo), ('augment_verify', report)]:
        state['steps'][step]['outputs'][str(path)] = RF.sha256(path)
    write_json(run / 'augmentation' / tag / 'state.json', state)
    # Change method-source provenance in tests without modifying shared production scripts.
    for module, name in [(RF.RS, 'redsea.py'), (RF.CT, 'codex.py')]:
        code = tmp_path / name
        code.write_text('method version 1')
        monkeypatch.setattr(module, '__file__', str(code))
    monkeypatch.setattr(RF, 'LOADED_CODE_SHA256', {
        'redsea_code': RF.sha256(Path(RF.RS.__file__)), 'typing_code': RF.sha256(Path(RF.CT.__file__)),
        'refresh_code': RF.sha256(Path(RF.__file__))})
    monkeypatch.setattr(RF.CT, 'cohort_rows', lambda: [{'image_name': tag}])
    return run, tag


def fake_output(run, tag, dest_dir):
    dest_dir.mkdir(parents=True, exist_ok=True)
    ids = pd.read_parquet(run / 'cells' / f'{tag}.parquet')['Object ID']
    frame = pd.DataFrame({'object_id': ids, 'qc_analysis_eligible': True,
                          'qc_spillover_context_eligible': True, 'qc_flags': ''})
    for comp in RF.RS.COMPARTMENTS:
        frame[f'{comp}__DAPI'] = 1.
        frame[f'{comp}__n_px'] = 10.
    frame.to_parquet(dest_dir / f'{tag}.parquet', index=False)
    record = {'image': tag, 'n_cells': len(ids), 'channels': ['DAPI'],
              'image_path': RF.CT.manifest(run)[tag]['image_path'], 'redsea': RF.CT.REDSEA,
              'spillover_keeps_no_nucleus': RF.RS.SPILLOVER_KEEPS_NO_NUCLEUS,
              'despeckle': RF.RS.DESPECKLE,
              'empty_area_null': {'method': RF.RS.EMPTY, 'per_channel': {'DAPI': {'null_p99': .5}}},
              'alignment_check': {'min_r': 1., 'per_channel': {'DAPI': {'r': 1.}}}}
    write_json(dest_dir / f'{tag}.json', record)
    return record


def install_processor(monkeypatch, mutate=None):
    calls = []
    def process(run, tag, force=False, dest_dir=None):
        assert force is True
        assert dest_dir != RF.CT.typing_dir(run) / 'redsea'
        calls.append(tag)
        result = fake_output(run, tag, dest_dir)
        if mutate:
            mutate(run, tag, dest_dir)
        return result
    monkeypatch.setattr(RF.RS, 'process', process)
    return calls


def update_checkpoint(run, tag, step, path):
    state_path = run / 'augmentation' / tag / 'state.json'
    state = json.loads(state_path.read_text())
    state['steps'][step]['outputs'][str(path)] = RF.sha256(path)
    if step == 'augment_detect':
        state['inputs']['outline_sha256'] = RF.sha256(path)
    write_json(state_path, state)


@pytest.mark.parametrize('unfinished', ['missing_report', 'missing_verify', 'missing_export', 'in_progress'])
def test_fullstate_required_before_submission(completed, monkeypatch, unfinished):
    run, tag = completed
    state_path = run / 'augmentation' / tag / 'state.json'
    state = json.loads(state_path.read_text())
    if unfinished == 'missing_report':
        (state_path.parent / 'verification.json').unlink()
    elif unfinished == 'in_progress':
        state['in_progress'] = {'id': 'augment_verify'}
    else:
        del state['steps']['augment_verify' if unfinished == 'missing_verify' else 'augment_J_context']
    write_json(state_path, state)
    calls = install_processor(monkeypatch)
    result = RF.refresh_ready(run, [tag], executor_factory=ThreadPoolExecutor)
    assert result == {'completed': {}, 'pending': [tag], 'failed': {}}
    assert calls == []


@pytest.mark.parametrize('key,rel', [('cells', 'cells/A.parquet'), ('geometry', 'cell_geojson/A.geojson'),
                                  ('tissue', 'inputs/tissue/A.geojson'),
                                  ('verification', 'augmentation/A/verification.json')])
def test_verified_input_hashes_are_enforced(completed, key, rel):
    run, tag = completed
    path = run / rel
    path.write_bytes(path.read_bytes() + b'changed')
    with pytest.raises(ValueError, match=f'verified {key} changed'):
        RF.ready_inputs(run, tag)


def test_copied_cache_forced_once_then_current_provenance_can_resume(completed, monkeypatch):
    run, tag = completed
    fake_output(run, tag, RF.CT.typing_dir(run) / 'redsea')
    calls = install_processor(monkeypatch)
    assert RF.refresh_one(run, tag)['status'] == 'refreshed'
    assert RF.refresh_one(run, tag)['status'] == 'current'
    assert calls == [tag]
    saved = json.loads(RF.provenance_path(run, tag).read_text())
    assert saved['inputs']['files']['cells']['sha256'] == RF.sha256(run / 'cells/A.parquet')
    assert saved['outputs']['parquet'] == RF.sha256(run / 'typing/redsea/A.parquet')


@pytest.mark.parametrize('changed', ['cells', 'geometry', 'tissue', 'redsea_code', 'image', 'cached_output'])
def test_stale_caches_recomputed(completed, monkeypatch, changed):
    run, tag = completed
    calls = install_processor(monkeypatch)
    RF.refresh_one(run, tag)
    if changed == 'cells':
        path = run / 'cells/A.parquet'
        frame = pd.read_parquet(path)
        frame.loc[0, 'Cell: DAPI: Mean'] = 6.
        frame.to_parquet(path, index=False)
        update_checkpoint(run, tag, 'augment_C_cells', path)
    elif changed in ('geometry', 'tissue'):
        path = run / ('cell_geojson/A.geojson' if changed == 'geometry' else 'inputs/tissue/A.geojson')
        value = json.loads(path.read_text())
        value['revision'] = 2
        write_json(path, value)
        update_checkpoint(run, tag, 'augment_C_geojson' if changed == 'geometry' else 'augment_detect', path)
    elif changed == 'redsea_code':
        Path(RF.RS.__file__).write_text('method version 2')
        # Simulate a restarted helper, which imports the new method before considering caches.
        RF.LOADED_CODE_SHA256['redsea_code'] = RF.sha256(Path(RF.RS.__file__))
    elif changed == 'image':
        Path(RF.CT.manifest(run)[tag]['image_path']).write_bytes(b'changed source image')
    else:
        path = run / 'typing/redsea/A.parquet'
        frame = pd.read_parquet(path)
        frame.loc[0, 'Cell__DAPI'] = 99.
        frame.to_parquet(path, index=False)
    assert RF.refresh_one(run, tag)['status'] == 'refreshed'
    assert calls == [tag, tag]


def test_input_change_during_processing_cannot_publish(completed, monkeypatch):
    run, tag = completed
    def mutate(run, tag, dest):
        tissue = run / 'inputs/tissue/A.geojson'
        write_json(tissue, {'features': [], 'new_revision': True})
        update_checkpoint(run, tag, 'augment_detect', tissue)
    install_processor(monkeypatch, mutate)
    with pytest.raises(ValueError, match='inputs changed during REDSEA'):
        RF.refresh_one(run, tag)
    assert not RF.provenance_path(run, tag).exists()
    assert not (run / 'typing/redsea/A.parquet').exists()


def test_live_code_change_requires_restart(completed, monkeypatch):
    run, tag = completed
    Path(RF.RS.__file__).write_text('source changed after its module was imported')
    calls = install_processor(monkeypatch)
    with pytest.raises(ValueError, match='changed since import; restart'):
        RF.refresh_one(run, tag)
    assert calls == []


@pytest.mark.parametrize('bad', ['duplicate', 'missing', 'extra', 'missing_column', 'bad_alignment'])
def test_invalid_redsea_outputs_not_published(completed, monkeypatch, bad):
    run, tag = completed
    def mutate(run, tag, dest):
        if bad == 'bad_alignment':
            path = dest / f'{tag}.json'
            rec = json.loads(path.read_text())
            rec['alignment_check']['per_channel']['DAPI']['r'] = .5
            write_json(path, rec)
            return
        path = dest / f'{tag}.parquet'
        frame = pd.read_parquet(path)
        if bad == 'duplicate':
            frame.loc[1, 'object_id'] = 'old'
        elif bad == 'missing':
            frame = frame.iloc[:1]
        elif bad == 'extra':
            frame.loc[1, 'object_id'] = 'unverified'
        else:
            frame = frame.drop(columns=['Cell__DAPI'])
        frame.to_parquet(path, index=False)
    install_processor(monkeypatch, mutate)
    with pytest.raises(ValueError):
        RF.refresh_one(run, tag)
    assert not RF.provenance_path(run, tag).exists()
    assert not (run / 'typing/redsea/A.parquet').exists()


def test_wait_schedules_later_verified_images_and_caps_parallelism(monkeypatch, tmp_path):
    tags = ['A', 'B', 'C', 'D', 'E', 'F']
    ready, done = {'A'}, set()
    lock = threading.Lock()
    active, maximum = 0, 0
    barrier = threading.Barrier(3)
    monkeypatch.setattr(RF, 'ready_inputs', lambda run, tag: {'tag': tag} if tag in ready else None)
    monkeypatch.setattr(RF, 'cache_current', lambda run, tag, snap: tag in done)
    def process(run, tag):
        nonlocal active, maximum
        with lock:
            assert tag in ready
            active += 1
            maximum = max(active, maximum)
        if tag == 'A':
            # A finishes before other augmentation checkpoints exist. The watcher must
            # discover newly verified images rather than exit or prequeue unfinished work.
            ready.update(tags[1:])
        if tag in ('B', 'C', 'D'):
            barrier.wait(timeout=3)
        with lock:
            done.add(tag)
            active -= 1
        return {'image': tag, 'status': 'refreshed'}
    monkeypatch.setattr(RF, 'refresh_one', process)
    result = RF.refresh_ready(tmp_path, tags, workers=3, watch=True, poll_seconds=.01,
                              executor_factory=ThreadPoolExecutor)
    assert maximum == 3
    assert result['pending'] == [] and result['failed'] == {}
    assert set(result['completed']) == set(tags)


def test_wait_with_nothing_ready_polls_until_verified(monkeypatch, tmp_path):
    ready = False
    def pause(seconds):
        nonlocal ready
        assert seconds == .01
        ready = True
    monkeypatch.setattr(RF.time, 'sleep', pause)
    monkeypatch.setattr(RF, 'ready_inputs', lambda run, tag: {'ready': True} if ready else None)
    monkeypatch.setattr(RF, 'refresh_one', lambda run, tag: {'status': 'refreshed'})
    monkeypatch.setattr(RF, 'cache_current', lambda *args: True)
    result = RF.refresh_ready(tmp_path, ['A'], watch=True, poll_seconds=.01,
                              executor_factory=ThreadPoolExecutor)
    assert result['completed'] == {'A': 'refreshed'}


def test_completed_input_change_while_watching_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(RF, 'ready_inputs', lambda *args: {'ready': True})
    monkeypatch.setattr(RF, 'refresh_one', lambda *args: {'status': 'refreshed'})
    monkeypatch.setattr(RF, 'cache_current', lambda *args: False)
    result = RF.refresh_ready(tmp_path, ['A'], executor_factory=ThreadPoolExecutor)
    assert 'changed after completion' in result['failed']['A']


def test_combined_nulls_wait_for_entire_fresh_cohort(completed, monkeypatch):
    run, tag = completed
    fake_output(run, tag, RF.CT.typing_dir(run) / 'redsea')
    dest = run / 'typing/empty_area_nulls.csv'
    dest.write_text('old nulls must not be republished')
    assert RF.publish_empty_nulls(run, [tag]) == [tag]
    assert dest.read_text() == 'old nulls must not be republished'
    install_processor(monkeypatch)
    RF.refresh_one(run, tag)
    assert RF.publish_empty_nulls(run, [tag]) == []
    rows = pd.read_csv(dest)
    assert rows[['image', 'marker']].values.tolist() == [[tag, 'DAPI']]
    assert (run / 'typing/redsea_refresh/empty_area_nulls.json').is_file()


def test_combining_nulls_rejects_concurrent_input_change(completed, monkeypatch):
    run, tag = completed
    install_processor(monkeypatch)
    RF.refresh_one(run, tag)
    real_write = RF.RS.write_empty_nulls
    def changing_write(directory):
        real_write(directory)
        Path(RF.CT.manifest(run)[tag]['image_path']).write_bytes(b'changed image')
    monkeypatch.setattr(RF.RS, 'write_empty_nulls', changing_write)
    with pytest.raises(ValueError, match='changed while combining'):
        RF.publish_empty_nulls(run, [tag])
    assert not (run / 'typing/empty_area_nulls.csv').exists()


def test_stale_combined_table_is_preserved_but_not_exposed(completed):
    run, tag = completed
    path = run / 'typing/empty_area_nulls.csv'
    path.parent.mkdir(parents=True)
    path.write_text('copied old nulls')
    archive = RF.retire_empty_nulls(run)
    assert not path.exists()
    assert (Path(archive) / path.name).read_text() == 'copied old nulls'
