#!/usr/bin/env python3
"""Refresh unchanged REDSEA measurements from a verified new-detection filter generation.

The raw augmentation reader stays frozen while its cohort jobs run. This adapter uses
the separate filter ancestry resolver and the existing REDSEA computation/output checks.
Scientific image-review release is an additional orchestration gate, not inferred here.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import refresh_ready_cells as RF
import new_detection_filter as NF

CT, RS = RF.CT, RF.RS
SCHEMA = RF.SCHEMA
sha256, signature, atomic_json = RF.sha256, RF.signature, RF.atomic_json
provenance_path, validate_outputs, cache_current = RF.provenance_path, RF.validate_outputs, RF.cache_current
retire_empty_nulls = RF.retire_empty_nulls
CODE_PATHS = {'redsea_code': Path(RS.__file__), 'typing_code': Path(CT.__file__),
              'refresh_code': Path(RF.__file__), 'filtered_refresh_code': Path(__file__)}
LOADED_CODE = {key: sha256(path) for key, path in CODE_PATHS.items()}


def ready_inputs(run: Path, tag: str) -> dict | None:
    snapshot = NF.accepted_inputs(run, tag)
    if snapshot is None:
        return None
    code = {key: signature(path) for key, path in CODE_PATHS.items()}
    for key, expected in LOADED_CODE.items():
        if code[key]['sha256'] != expected:
            raise ValueError(f'{tag}: {key} changed since import; restart the filtered refresh')
    return {**snapshot, 'files': {**snapshot['files'], **code}}


def refresh_one(run: Path, tag: str) -> dict:
    before = ready_inputs(run, tag)
    if before is None:
        raise ValueError(f'{tag}: filtered generation is incomplete')
    if cache_current(run, tag, before):
        return {'image': tag, 'status': 'current'}
    dest = CT.typing_dir(run) / 'redsea'
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f'.filtered-redsea-{tag}-', dir=dest.parent) as folder:
        staging = Path(folder)
        RS.process(run, tag, force=True, dest_dir=staging)
        outputs = validate_outputs(run, tag, staging, before)
        if ready_inputs(run, tag) != before:
            raise ValueError(f'{tag}: filtered inputs changed during REDSEA; publication withheld')
        dest.mkdir(parents=True, exist_ok=True)
        for suffix in ('parquet', 'json'):
            os.replace(staging / f'{tag}.{suffix}', dest / f'{tag}.{suffix}')
        if ready_inputs(run, tag) != before:
            raise ValueError(f'{tag}: filtered inputs changed during publication; cache remains invalid')
        atomic_json(provenance_path(run, tag), {'schema': SCHEMA, 'inputs': before, 'outputs': outputs,
                    'completed_at_utc': datetime.now(timezone.utc).isoformat()})
    return {'image': tag, 'status': 'refreshed'}


def refresh_ready(run: Path, tags: list[str], *, workers: int = 3) -> dict:
    if workers not in (1, 2, 3):
        raise ValueError('workers must be 1..3')
    pending, completed, failed = [], {}, {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        active = {}
        for tag in dict.fromkeys(tags):
            try:
                if ready_inputs(run, tag) is None:
                    pending.append(tag)
                else:
                    active[pool.submit(refresh_one, run, tag)] = tag
            except Exception as exc:
                failed[tag] = str(exc)
        for future in as_completed(active):
            tag = active[future]
            try:
                completed[tag] = future.result()['status']
                print(f'{tag}: {completed[tag]}', flush=True)
            except BaseException as exc:
                failed[tag] = str(exc)
        for tag in completed:
            try:
                snapshot = ready_inputs(run, tag)
                if snapshot is None or not cache_current(run, tag, snapshot):
                    raise ValueError('Filtered inputs or REDSEA outputs changed after completion')
            except Exception as exc:
                failed[tag] = str(exc)
    return {'completed': completed, 'pending': pending, 'failed': failed}


def publish_empty_nulls(run: Path, cohort: list[str]) -> list[str]:
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
    with tempfile.TemporaryDirectory(prefix='.filtered-redsea-nulls-', dir=typing) as folder:
        stage = Path(folder) / 'redsea'
        stage.mkdir()
        for tag in cohort:
            raw = (typing / 'redsea' / f'{tag}.json').read_bytes()
            records[tag] = hashlib.sha256(raw).hexdigest()
            (stage / f'{tag}.json').write_bytes(raw)
        RS.write_empty_nulls(stage)
        nulls = stage.parent / 'empty_area_nulls.csv'
        if not nulls.is_file():
            raise ValueError('No empty-slide null measurements for the filtered cohort')
        for tag in cohort:
            if ready_inputs(run, tag) != inputs[tag] or not cache_current(run, tag, inputs[tag]):
                raise ValueError(f'{tag}: inputs or REDSEA outputs changed during null-table assembly')
            if sha256(typing / 'redsea' / f'{tag}.json') != records[tag]:
                raise ValueError(f'{tag}: REDSEA record changed during null-table assembly')
        os.replace(nulls, typing / 'empty_area_nulls.csv')
    atomic_json(typing / 'redsea_refresh' / 'empty_area_nulls.json',
                {'schema': SCHEMA, 'redsea_record_sha256': records,
                 'csv_sha256': sha256(typing / 'empty_area_nulls.csv')})
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--images', required=True, nargs='+')
    parser.add_argument('--workers', choices=(1, 2, 3), type=int, default=3)
    args = parser.parse_args()
    run = args.run.resolve()
    cohort = [r['image_name'] for r in CT.cohort_rows()]
    tags = cohort if args.images == ['all'] else args.images
    if set(tags) - set(cohort) or len(tags) != len(set(tags)):
        parser.error('images must be unique cohort names, or all')
    if not (run / 'filter/generation.json').is_file():
        parser.error('--run must contain a verified filter generation')
    lock_path = CT.typing_dir(run) / 'redsea_refresh/driver.lock'
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = refresh_ready(run, tags, workers=args.workers)
        if not result['failed']:
            result['empty_nulls_pending'] = publish_empty_nulls(run, cohort)
        print(json.dumps(result, indent=2), flush=True)
        return int(bool(result['failed']))


if __name__ == '__main__':
    raise SystemExit(main())
