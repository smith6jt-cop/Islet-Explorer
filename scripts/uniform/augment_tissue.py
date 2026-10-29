#!/usr/bin/env python3
"""Augment a COPY of a completed uniform run without rebuilding its cells or endocrine objects.

    python scripts/uniform/augment_tissue.py --source-run <baseline> --out <copied-run> \\
        --outlines <reviewed-content-outlines> --images ND_0112 --stages detect measure export verify

The copy must have independent qpdata files. This driver never invokes tissue_import, the full
detection step, islet fitting, step 1/1b, or islet naming. Its own checkpoints are kept under
augmentation/<image>/ and provenance/<image>/augment_*.json. Repeating the same command resumes
completed stages. Changed inputs or externally modified qpdata are refused, never overwritten.
Use --dry-run to inspect stage order and inputs without creating files or running QuPath.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import threading

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_uniform as RU

STAGES = ('detect', 'measure', 'export', 'verify')
EXPORT_STEPS = ('3', '8', 'C_cells', 'C_geojson', 'J_cells', 'J_context', 'J_annotations', 'J_islets')
QUANTITATIVE_PREFIXES = ('Cell:', 'Nucleus:', 'Cytoplasm:', 'Membrane:')
GPU_STAGE_LOCK = threading.Lock()
DETECTION_CONTEXT_UM = 30.0


def read_json(path: Path):
    return json.loads(path.read_text())


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + '.part')
    part.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(part, path)


def digest_json(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def qpdata(run: Path, tag: str) -> Path:
    entries = read_json(run / 'qupath' / 'project.qpproj')['images']
    matches = [e for e in entries if e['imageName'] == tag]
    if len(matches) != 1:
        raise ValueError(f'{tag}: expected exactly one project entry')
    return run / 'qupath' / 'data' / str(matches[0]['entryID']) / 'data.qpdata'


def check_independent(source: Path, target: Path) -> None:
    if not source.is_file() or not target.is_file():
        raise ValueError(f'Missing source or copied qpdata: {source}, {target}')
    if source.resolve() == target.resolve() or os.path.samefile(source, target):
        raise ValueError('Output qpdata must be an independent copy, not a symlink or hardlink to source')


def stage_steps(stage: str) -> list[RU.Step]:
    cfg = '{OUT}/augmentation/{TAG}/config.json'
    if stage == 'detect':
        return [RU.Step('augment_detect', 'qupath', 'qupath/augment_tissue.groovy', (cfg,),
                        'footer-saves', 96, 12, 240, 'v3_result', gpu=True)]
    if stage == 'measure':
        return [RU.Step('augment_measure', 'qupath', 'qupath/measure_augmented_cells.groovy', (cfg,),
                        'footer-saves', 128, 24, 240, 'v3_result')]
    if stage == 'verify':
        return [RU.Step('augment_verify', 'qupath', 'qupath/verify_tissue_augmentation.groovy', (cfg,),
                        'never', 96, 8, 90, 'v3_result')]
    if stage != 'export':
        raise ValueError(stage)
    result = []
    for name in EXPORT_STEPS:
        step = replace(RU.BY_ID[name], id='augment_' + name, after=(), post='')
        if name in ('3', '8'):
            step = replace(step, script='qupath/refresh_augmented_context.groovy', args=(cfg, name))
        result.append(step)
    return result


def post_export(runner: RU.Runner, tag: str, step: RU.Step) -> list[Path]:
    """Convert exports using the existing table contract without touching baseline provenance."""
    import pyarrow.csv as pcsv
    import pyarrow.parquet as pq
    name = step.id.removeprefix('augment_')
    original = RU.BY_ID.get(name)
    rec = read_json(runner.record(tag, step.id))
    if original is None or not original.post:
        return [Path(p['path']) for p in rec.get('outputs', []) if p.get('exists')]
    tsv = Path(runner.subst(original.args[1], tag))
    dest = runner.out / RU.POST_TARGET[(original.id, original.post)].replace('{TAG}', tag)
    table = pcsv.read_csv(tsv, parse_options=pcsv.ParseOptions(delimiter='\t'),
                          convert_options=pcsv.ConvertOptions(strings_can_be_null=True))
    ids = table.column('Object ID').to_pandas()
    if table.num_rows != rec['counts']['detections'] or ids.isna().any() or ids.duplicated().any():
        raise ValueError(f'{tag}/{step.id}: incomplete or duplicated export IDs')
    if original.post == 'cells_parquet' and 'Cell: DAPI: Mean' not in table.column_names:
        raise ValueError(f'{tag}/{step.id}: missing DAPI measurements')
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + '.part')
    pq.write_table(table, part)
    os.replace(part, dest)
    # Keep the TSV until the whole stage is checkpointed, so crash recovery can replay conversion.
    return [dest]


def verify_exports(source: Path, out: Path, tag: str, expected_new: set[str]) -> dict:
    """Compare every original raw quantitative value, IDs and serialized cell/endocrine geometry."""
    import numpy as np
    import pandas as pd
    old = pd.read_parquet(source / 'cells' / f'{tag}.parquet').set_index('Object ID')
    new = pd.read_parquet(out / 'cells' / f'{tag}.parquet').set_index('Object ID')
    if not old.index.is_unique or not new.index.is_unique:
        raise ValueError(f'{tag}: duplicate export IDs')
    if set(new.index) != set(old.index) | expected_new or set(old.index) & expected_new:
        raise ValueError(f'{tag}: exported IDs do not equal original plus accepted additions')
    columns = [c for c in old.columns if c.startswith(QUANTITATIVE_PREFIXES) or c.startswith('Centroid ')]
    if set(columns) - set(new.columns):
        raise ValueError(f'{tag}: original measurement columns disappeared')
    left, right = old[columns], new.loc[old.index, columns]
    for c in columns:
        a, b = left[c].to_numpy(), right[c].to_numpy()
        if not np.array_equal(a, b, equal_nan=True):
            raise ValueError(f'{tag}: original measurement changed: {c}')
    def features(path):
        raw = read_json(path)
        return raw['features'] if isinstance(raw, dict) else raw
    old_geo = {f['id']: f for f in features(source / 'cell_geojson' / f'{tag}.geojson')}
    new_geo = {f['id']: f for f in features(out / 'cell_geojson' / f'{tag}.geojson')}
    if set(new_geo) != set(new.index):
        raise ValueError(f'{tag}: GeoJSON and measurements have different IDs')
    for oid, f in old_geo.items():
        for key in ('geometry', 'nucleusGeometry'):
            if f.get(key) != new_geo[oid].get(key):
                raise ValueError(f'{tag}: original {key} changed for {oid}')
    def endocrine(path):
        return {f['id']: f for f in features(path)
                if (f.get('properties', {}).get('classification') or {}).get('name')
                in ('Islet', 'EndocrineCluster', 'EndocrineSingle')}
    old_endo = endocrine(source / 'islets' / f'{tag}.geojson')
    new_endo = endocrine(out / 'islets' / f'{tag}.geojson')
    if old_endo != new_endo:
        raise ValueError(f'{tag}: accepted endocrine objects changed')
    for folder in ('final_cells', 'context'):
        table = pd.read_parquet(out / folder / f'{tag}.parquet', columns=['Object ID'])
        if table['Object ID'].duplicated().any() or set(table['Object ID']) != set(new.index):
            raise ValueError(f'{tag}: {folder} inventory differs')
    return {'original_cells': len(old), 'new_cells': len(expected_new), 'total_cells': len(new),
            'raw_columns_preserved': len(columns), 'endocrine_objects_preserved': len(old_endo)}


def process(runner: RU.Runner, source: Path, outlines: Path, tag: str, stages: list[str],
            detection_context_um: float = DETECTION_CONTEXT_UM) -> dict:
    if not math.isfinite(detection_context_um) or detection_context_um <= 0:
        raise ValueError('detection_context_um must be finite and positive')
    area = runner.out / 'augmentation' / tag
    state_path, config_path = area / 'state.json', area / 'config.json'
    src_qp, dst_qp = qpdata(source, tag), qpdata(runner.out, tag)
    check_independent(src_qp, dst_qp)
    outline = outlines / f'{tag}.geojson'
    signature = {'source_run': str(source), 'source_qpdata_sha256': RU.sha256(src_qp),
                 'outline': str(outline), 'outline_sha256': RU.sha256(outline),
                 'detection_context_um': float(detection_context_um)}
    current = RU.sha256(dst_qp)
    if state_path.exists():
        state = read_json(state_path)
        if state['inputs'] != signature:
            raise ValueError(f'{tag}: augmentation inputs changed; use a fresh copy for a changed proposal')
    else:
        if current != signature['source_qpdata_sha256']:
            raise ValueError(f'{tag}: initial output qpdata does not match source')
        state = {'schema': 'uniform.tissue_augmentation/1', 'inputs': signature,
                 'last_qpdata_sha256': current, 'steps': {}}
        atomic_json(state_path, state)
    config = {'outline': str(outline), 'inventory': str(area / 'original_inventory.json'),
              'new_ids': str(area / 'new_object_ids.json'),
              'detection_context_um': float(detection_context_um),
              'invariants_script': str(HERE / 'qupath' / 'augmentation_invariants.groovy')}
    atomic_json(config_path, config)
    all_steps = [step for stage in STAGES for step in stage_steps(stage)]
    requested = {step.id for stage in stages for step in stage_steps(stage)}
    for i, step in enumerate(all_steps):
        if step.id not in requested:
            continue
        needed = [s.id for s in all_steps[:i]]
        missing = [s for s in needed if s not in state['steps']]
        if missing:
            raise ValueError(f'{tag}/{step.id}: missing augmentation stages: {missing}')
        script_hash = RU.sha256(HERE / step.script)
        helper_hash = RU.sha256(Path(config['invariants_script']))
        source_name = step.id.removeprefix('augment_')
        delegated_hash = RU.sha256(HERE / RU.BY_ID[source_name].script) if source_name in ('3', '8') else None
        fingerprint = digest_json({'inputs': signature, 'script': script_hash, 'helper': helper_hash, 'delegated': delegated_hash,
                                   'footer': RU.sha256(RU.FOOTER)})
        if step.id in state['steps']:
            saved = state['steps'][step.id]
            if saved['fingerprint'] != fingerprint:
                raise ValueError(f'{tag}/{step.id}: code changed since checkpoint; inspect before resuming')
            for path, sha in saved['outputs'].items():
                if not Path(path).is_file() or RU.sha256(Path(path)) != sha:
                    raise ValueError(f'{tag}/{step.id}: checkpoint output changed: {path}')
            print(f'{tag}/{step.id}: already complete', flush=True)
            continue
        pending = state.get('in_progress')
        record_path = runner.record(tag, step.id)
        recover = False
        if pending and pending['id'] == step.id and pending['fingerprint'] == fingerprint and record_path.is_file():
            rec = read_json(record_path)
            after = rec.get('qpdata', {}).get('after', {})
            observed = RU.stat_of(dst_qp)
            recover = (rec.get('status') == 'pass' and rec.get('script', {}).get('source_sha256') == script_hash
                       and observed.get('size') == after.get('size') and observed.get('mtime_ns') == after.get('mtime_ns'))
        if not recover:
            if RU.sha256(dst_qp) != state['last_qpdata_sha256']:
                raise ValueError(f'{tag}: qpdata changed outside a completed checkpoint; inspect before recovery')
            if record_path.exists():
                runner.retire(record_path, 'augmentation-retry')
            state['in_progress'] = {'id': step.id, 'fingerprint': fingerprint}
            atomic_json(state_path, state)
            config['deadline_epoch_ms'] = int((time.time() + step.timeout_min * 60) * 1000)
            atomic_json(config_path, config)
            if step.gpu:
                with GPU_STAGE_LOCK:
                    ok = runner.run_qupath(tag, step)
            else:
                ok = runner.run_qupath(tag, step)
            if not ok:
                raise RuntimeError(f'{tag}/{step.id}: failed; see {runner.out}/logs/{tag}/{step.id}.log')
        else:
            print(f'{tag}/{step.id}: recovering saved successful stage', flush=True)
        outputs = post_export(runner, tag, step)
        if step.id == 'augment_detect':
            dest_outline = runner.out / 'inputs' / 'tissue' / f'{tag}.geojson'
            dest_outline.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(outline, dest_outline)
            outputs.append(dest_outline)
        if step.id == 'augment_verify':
            report = verify_exports(source, runner.out, tag, set(read_json(Path(config['new_ids']))))
            report_path = area / 'verification.json'
            atomic_json(report_path, report)
            outputs.append(report_path)
        state['last_qpdata_sha256'] = RU.sha256(dst_qp)
        state['steps'][step.id] = {'fingerprint': fingerprint,
                                  'outputs': {str(p): RU.sha256(p) for p in outputs}}
        state.pop('in_progress', None)
        atomic_json(state_path, state)
        print(f'{tag}/{step.id}: complete', flush=True)
    if RU.sha256(dst_qp) != state['last_qpdata_sha256']:
        raise ValueError(f'{tag}: qpdata changed after the final checkpoint')
    return state


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--source-run', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--outlines', required=True, type=Path)
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--stages', nargs='+', choices=STAGES, default=list(STAGES))
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--detection-context-um', type=float, default=DETECTION_CONTEXT_UM,
                    help='technical detection-ROI context buffer, clipped only to image bounds (default: 30 um); '
                         'acceptance still uses the unbuffered added-Tissue centroid and full final-Tissue containment')
    ap.add_argument('--workers', type=int, choices=(1, 2, 3), default=1,
                    help='independent images in parallel; GPU stages remain serialized; at most three 128-GB JVMs')
    a = ap.parse_args()
    if not math.isfinite(a.detection_context_um) or a.detection_context_um <= 0:
        ap.error('--detection-context-um must be finite and positive')
    source, out, outlines = a.source_run.resolve(), a.out.resolve(), a.outlines.resolve()
    if source == out:
        ap.error('--out must be a separate copy of --source-run')
    if list(dict.fromkeys(a.stages)) != a.stages or a.stages != sorted(a.stages, key=STAGES.index):
        ap.error('--stages must be unique and in detect measure export verify order')
    cohort_path = HERE / 'cohort.csv'
    rows = list(csv.DictReader(cohort_path.open()))
    tags = [r['image_name'] for r in rows] if a.images == ['all'] else a.images
    if set(tags) - {r['image_name'] for r in rows}:
        ap.error('unknown images')
    runner = RU.Runner(out, rows, cohort_path)
    if a.dry_run:
        for tag in tags:
            check_independent(qpdata(source, tag), qpdata(out, tag))
            if not (outlines / f'{tag}.geojson').is_file():
                raise ValueError(f'{tag}: missing outline')
        print(json.dumps({'source_run': str(source), 'out': str(out), 'outlines': str(outlines), 'images': tags,
                          'detection_context_um': a.detection_context_um,
                          'steps': [s.id for stage in a.stages for s in stage_steps(stage)]}, indent=2))
        return 0
    lock_path = out / 'augmentation' / 'driver.lock'
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with ThreadPoolExecutor(a.workers) as pool:
            list(pool.map(lambda tag: process(runner, source, outlines, tag, a.stages, a.detection_context_um), tags))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
