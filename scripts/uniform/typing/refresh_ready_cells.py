#!/usr/bin/env python3
"""Refresh REDSEA as each image finishes the full tissue-augmentation checkpoint chain.

    python scripts/uniform/typing/refresh_ready_cells.py --run <candidate> --images all --wait

At most three image processes overlap the augmentation driver. Copied REDSEA files are never
accepted without this helper's current per-image provenance. Outputs are calculated in scratch
space and published only after their inventory and unchanged, verified inputs have been checked.
Without --wait, drain currently ready work and report images still awaiting augmentation.
The cohort-wide empty-slide null table is published only when the whole cohort is current.
"""
from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT
import redsea_cells as RS

SCHEMA = 'uniform.redsea_refresh/1'
# augment_tissue.py's full detect, measure, export, verify chain. A pilot/detection-only
# state, an uncheckpointed verify, or an in-progress stage must never release an image.
REQUIRED_STEPS = ('augment_detect', 'augment_measure', 'augment_3', 'augment_8',
                  'augment_C_cells', 'augment_C_geojson', 'augment_J_cells',
                  'augment_J_context', 'augment_J_annotations', 'augment_J_islets', 'augment_verify')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def file_stat(path: Path) -> dict:
    st = path.stat()
    return {'path': str(path.resolve()), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns,
            'ctime_ns': st.st_ctime_ns, 'inode': st.st_ino}


def signature(path: Path) -> dict:
    before = file_stat(path)
    result = {**before, 'sha256': sha256(path)}
    if file_stat(path) != before:
        raise ValueError(f'Input changed while being hashed: {path}')
    return result


# A long-lived watcher/process pool keeps imported Python code in memory. If a source
# changes on disk, require a restart rather than crediting old code with a new hash.
LOADED_CODE_SHA256 = {key: sha256(path) for key, path in {
    'redsea_code': Path(RS.__file__), 'typing_code': Path(CT.__file__), 'refresh_code': Path(__file__)
}.items()}


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix=path.name + '.', delete=False) as f:
        temporary = Path(f.name)
        json.dump(value, f, indent=2)
        f.write('\n')
    os.replace(temporary, path)


def provenance_path(run: Path, tag: str) -> Path:
    return CT.typing_dir(run) / 'redsea_refresh' / f'{tag}.json'


def ready_inputs(run: Path, tag: str) -> dict | None:
    """Return the verified input snapshot, None for incomplete augmentation, or fail stale checkpoints."""
    area = run / 'augmentation' / tag
    state_path, report_path = area / 'state.json', area / 'verification.json'
    if not state_path.is_file() or not report_path.is_file():
        return None
    state_bytes = state_path.read_bytes()
    state = json.loads(state_bytes)
    if state.get('schema') != 'uniform.tissue_augmentation/1':
        raise ValueError(f'{tag}: unsupported augmentation state schema')
    steps = state.get('steps', {})
    if state.get('in_progress') or set(REQUIRED_STEPS) - set(steps):
        return None
    for name in REQUIRED_STEPS:
        if not steps[name].get('fingerprint') or not isinstance(steps[name].get('outputs'), dict):
            raise ValueError(f'{tag}: malformed completed checkpoint {name}')
    needed = {
        'cells': ('augment_C_cells', run / 'cells' / f'{tag}.parquet'),
        'geometry': ('augment_C_geojson', run / 'cell_geojson' / f'{tag}.geojson'),
        'tissue': ('augment_detect', run / 'inputs' / 'tissue' / f'{tag}.geojson'),
        'verification': ('augment_verify', report_path),
    }
    files = {}
    for key, (step, path) in needed.items():
        expected = [v for p, v in steps[step]['outputs'].items() if Path(p).resolve() == path.resolve()]
        if len(expected) != 1:
            raise ValueError(f'{tag}: {step} does not checkpoint {path}')
        files[key] = signature(path)
        if files[key]['sha256'] != expected[0]:
            raise ValueError(f'{tag}: verified {key} changed after augmentation checkpoint')
    report = json.loads(report_path.read_text())
    total = report.get('total_cells')
    if not isinstance(total, int) or total < 0 or total != report.get('original_cells', -1) + report.get('new_cells', -1):
        raise ValueError(f'{tag}: invalid augmentation verification counts')
    if files['tissue']['sha256'] != state.get('inputs', {}).get('outline_sha256'):
        raise ValueError(f'{tag}: Tissue does not match the verified augmentation outline')
    entry = CT.manifest(run)[tag]
    server = run / 'qupath' / 'data' / str(entry['entry_id']) / 'server.json'
    files['server'] = signature(server)
    files['redsea_code'] = signature(Path(RS.__file__))
    files['typing_code'] = signature(Path(CT.__file__))
    files['refresh_code'] = signature(Path(__file__))
    for key, expected in LOADED_CODE_SHA256.items():
        if files[key]['sha256'] != expected:
            raise ValueError(f'{tag}: {key} changed since import; restart the watcher before processing')
    names = CT.channels(run, entry['entry_id'])
    snapshot = {'run': str(run.resolve()), 'image': tag, 'files': files, 'entry': entry,
                # TIFFs are immutable source images; record their identity/size/timestamps rather
                # than rehashing tens of GB every time an augmentation image becomes ready.
                'source_image': file_stat(Path(entry['image_path'])), 'channels': names,
                'augmentation_state_sha256': hashlib.sha256(state_bytes).hexdigest(), 'total_cells': total}
    if state_path.read_bytes() != state_bytes or sha256(report_path) != files['verification']['sha256']:
        raise ValueError(f'{tag}: augmentation checkpoint changed during input validation')
    return snapshot


def validate_outputs(run: Path, tag: str, directory: Path, inputs: dict) -> dict:
    """Check that the complete REDSEA inventory and method match this verified image."""
    ids = pd.read_parquet(run / 'cells' / f'{tag}.parquet', columns=['Object ID'])['Object ID']
    frame = pd.read_parquet(directory / f'{tag}.parquet')
    names = inputs['channels']
    required = {'object_id', 'qc_analysis_eligible', 'qc_spillover_context_eligible', 'qc_flags'}
    required.update(f'{comp}__{c}' for comp in RS.COMPARTMENTS for c in [*names, 'n_px'])
    if required - set(frame.columns):
        raise ValueError(f'{tag}: incomplete REDSEA columns: {sorted(required - set(frame.columns))}')
    actual = frame['object_id']
    if ids.isna().any() or ids.duplicated().any() or actual.isna().any() or actual.duplicated().any():
        raise ValueError(f'{tag}: null or duplicate input/REDSEA IDs')
    if len(ids) != inputs['total_cells'] or set(actual) != set(ids):
        raise ValueError(f'{tag}: REDSEA IDs differ from verified cells')
    record = json.loads((directory / f'{tag}.json').read_text())
    if (record.get('image') != tag or record.get('n_cells') != len(ids) or record.get('channels') != names
            or record.get('image_path') != inputs['entry']['image_path'] or record.get('redsea') != CT.REDSEA
            or record.get('spillover_keeps_no_nucleus') != RS.SPILLOVER_KEEPS_NO_NUCLEUS
            or any(record.get('despeckle', {}).get(k) != v for k, v in RS.DESPECKLE.items())
            or record.get('empty_area_null', {}).get('method') != RS.EMPTY):
        raise ValueError(f'{tag}: REDSEA metadata does not match current input/method')
    alignment = record.get('alignment_check', {})
    per_channel = alignment.get('per_channel', {})
    if (set(per_channel) != set(names) or not alignment.get('min_r', float('nan')) >= RS.ALIGN_MIN_R
            or any(not v.get('r', float('nan')) >= RS.ALIGN_MIN_R for v in per_channel.values())):
        raise ValueError(f'{tag}: REDSEA alignment validation did not pass')
    nulls = record['empty_area_null'].get('per_channel')
    if not isinstance(nulls, dict) or (nulls and set(nulls) != set(names)):
        raise ValueError(f'{tag}: incomplete empty-slide null record')
    return {name: sha256(directory / f'{tag}.{name}') for name in ('parquet', 'json')}


def cache_current(run: Path, tag: str, inputs: dict) -> bool:
    path = provenance_path(run, tag)
    if not path.is_file():
        return False
    try:
        saved = json.loads(path.read_text())
        if saved.get('schema') != SCHEMA or saved.get('inputs') != inputs:
            return False
        return saved.get('outputs') == validate_outputs(run, tag, CT.typing_dir(run) / 'redsea', inputs)
    except (OSError, ValueError, KeyError):
        return False


def refresh_one(run: Path, tag: str) -> dict:
    before = ready_inputs(run, tag)
    if before is None:
        raise ValueError(f'{tag}: full augmentation verification is no longer ready')
    if cache_current(run, tag, before):
        return {'image': tag, 'status': 'current'}
    dest = CT.typing_dir(run) / 'redsea'
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f'.redsea-refresh-{tag}-', dir=dest.parent) as folder:
        staging = Path(folder)
        RS.process(run, tag, force=True, dest_dir=staging)
        outputs = validate_outputs(run, tag, staging, before)
        if ready_inputs(run, tag) != before:
            raise ValueError(f'{tag}: inputs changed during REDSEA; results not published')
        dest.mkdir(parents=True, exist_ok=True)
        for suffix in ('parquet', 'json'):
            os.replace(staging / f'{tag}.{suffix}', dest / f'{tag}.{suffix}')
        if ready_inputs(run, tag) != before:
            raise ValueError(f'{tag}: inputs changed while publishing REDSEA; cache remains invalid')
        atomic_json(provenance_path(run, tag), {'schema': SCHEMA, 'inputs': before, 'outputs': outputs,
                    'completed_at_utc': datetime.now(timezone.utc).isoformat()})
    return {'image': tag, 'status': 'refreshed'}


def publish_empty_nulls(run: Path, cohort: list[str]) -> list[str]:
    """Never let RS.write_empty_nulls combine fresh records with copied, stale donor caches."""
    inputs, pending = {}, []
    for tag in cohort:
        snap = ready_inputs(run, tag)
        if snap is None or not cache_current(run, tag, snap):
            pending.append(tag)
        else:
            inputs[tag] = snap
    if pending:
        return pending
    typing = CT.typing_dir(run)
    records = {}
    with tempfile.TemporaryDirectory(prefix='.redsea-nulls-', dir=typing) as folder:
        stage = Path(folder) / 'redsea'
        stage.mkdir()
        for tag in cohort:
            raw = (typing / 'redsea' / f'{tag}.json').read_bytes()
            records[tag] = hashlib.sha256(raw).hexdigest()
            (stage / f'{tag}.json').write_bytes(raw)
        RS.write_empty_nulls(stage)
        nulls = stage.parent / 'empty_area_nulls.csv'
        if not nulls.is_file():
            raise ValueError('No empty-slide null measurements were produced for the verified cohort')
        for tag in cohort:
            if ready_inputs(run, tag) != inputs[tag] or not cache_current(run, tag, inputs[tag]):
                raise ValueError(f'{tag}: inputs or REDSEA outputs changed while combining empty-slide nulls')
            if sha256(typing / 'redsea' / f'{tag}.json') != records[tag]:
                raise ValueError(f'{tag}: REDSEA record changed while combining empty-slide nulls')
        os.replace(nulls, typing / 'empty_area_nulls.csv')
    atomic_json(typing / 'redsea_refresh' / 'empty_area_nulls.json',
                {'schema': SCHEMA, 'redsea_record_sha256': records,
                 'csv_sha256': sha256(typing / 'empty_area_nulls.csv')})
    return []


def retire_empty_nulls(run: Path) -> str | None:
    """Keep a copied/prior combined table out of downstream reach while the refresh is incomplete."""
    typing = CT.typing_dir(run)
    paths = [typing / 'empty_area_nulls.csv', typing / 'redsea_refresh' / 'empty_area_nulls.json']
    existing = [p for p in paths if p.exists()]
    if not existing:
        return None
    archive = typing / 'redsea_refresh' / 'superseded'
    archive.mkdir(parents=True, exist_ok=True)
    folder = Path(tempfile.mkdtemp(prefix='empty-nulls-', dir=archive))
    for path in existing:
        os.replace(path, folder / path.name)
    return str(folder)


def refresh_ready(run: Path, tags: list[str], *, workers: int = 3, watch: bool = False,
                  poll_seconds: float = 15, executor_factory=ProcessPoolExecutor) -> dict:
    """Submit only verified images, with at most workers futures in flight (no speculative queue)."""
    if workers not in (1, 2, 3) or not 0 < poll_seconds <= 60:
        raise ValueError('workers must be 1..3 and poll_seconds must be in (0, 60]')
    pending, done, failed = list(dict.fromkeys(tags)), {}, {}
    active = {}
    with executor_factory(max_workers=workers) as pool:
        while pending or active:
            for tag in list(pending):
                if len(active) >= workers or failed:
                    break
                try:
                    snap = ready_inputs(run, tag)
                except Exception as exc:
                    failed[tag] = str(exc)
                    break
                if snap is not None:
                    active[pool.submit(refresh_one, run, tag)] = tag
                    pending.remove(tag)
            if not active:
                if failed or not watch or not pending:
                    break
                print(f'Awaiting full augmentation verification: {", ".join(pending)}', flush=True)
                time.sleep(poll_seconds)
                continue
            finished, _ = wait(active, timeout=poll_seconds, return_when=FIRST_COMPLETED)
            for future in finished:
                tag = active.pop(future)
                try:
                    done[tag] = future.result()['status']
                    print(f'{tag}: {done[tag]}', flush=True)
                except BaseException as exc:  # RS can raise SystemExit; report the failed donor.
                    failed[tag] = str(exc)
                    print(f'{tag}: FAILED: {exc}', flush=True)
        # Images completed early in --wait must still be current when the watcher exits.
        for tag in done:
            try:
                snap = ready_inputs(run, tag)
                if snap is None or not cache_current(run, tag, snap):
                    raise ValueError('verified inputs or REDSEA outputs changed after completion')
            except Exception as exc:
                failed[tag] = str(exc)
    return {'completed': done, 'pending': pending, 'failed': failed}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', required=True, type=Path)
    ap.add_argument('--images', required=True, nargs='+')
    ap.add_argument('--workers', type=int, choices=(1, 2, 3), default=3)
    ap.add_argument('--wait', action='store_true', help='watch until every requested image is verified and refreshed')
    ap.add_argument('--poll-seconds', type=float, default=15)
    a = ap.parse_args()
    if not 0 < a.poll_seconds <= 60:
        ap.error('--poll-seconds must be in (0, 60]')
    cohort = [r['image_name'] for r in CT.cohort_rows()]
    tags = cohort if a.images == ['all'] else a.images
    if set(tags) - set(cohort) or len(tags) != len(set(tags)):
        ap.error('images must be unique cohort names, or all')
    run = a.run.resolve()
    if not (run / 'project_manifest.json').is_file():
        ap.error('--run must name an existing uniform run')
    lock_path = CT.typing_dir(run) / 'redsea_refresh' / 'driver.lock'
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        retired = retire_empty_nulls(run)
        result = refresh_ready(run, tags, workers=a.workers, watch=a.wait, poll_seconds=a.poll_seconds)
        result['previous_empty_nulls_archive'] = retired
        if not result['failed']:
            result['empty_nulls_pending'] = publish_empty_nulls(run, cohort)
        print(json.dumps(result, indent=2), flush=True)
        return 1 if result['failed'] else 0


if __name__ == '__main__':
    sys.exit(main())
