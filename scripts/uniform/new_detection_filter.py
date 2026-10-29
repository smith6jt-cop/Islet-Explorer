#!/usr/bin/env python3
"""Create an independent, verified generation using explicit reviewed NEW-cell decisions.

No scientific criterion is inferred here. Never changes the raw augmentation run.
Use --seal-only to obtain the source_seal required in each decisions image record.
Decisions schema uniform.new_detection_decisions/1:
  rule: {id, version, parameters, decision_basis}
  images: {TAG: {source_seal: <source_seal()>, review: {status: "reviewed", reviewer,
    evidence: [{path, sha256}]}, decisions: [{object_id, decision: "accept"|"reject"|"defer", reason}]}}

Completion is filter/TAG/manifest.json, written last. accepted_inputs() validates it and
its immutable ancestry; it does not itself grant scientific/visual-review release.
Established structure reuse is separately certified by --certify-structures.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / 'typing'))
import augment_tissue as AU
import codex_typing as CT

SCHEMA = 'uniform.new_detection_filter/1'
DECISIONS_SCHEMA = 'uniform.new_detection_decisions/1'
SCRIPT = HERE / 'qupath/filter_new_detections.groovy'
HELPER = HERE / 'qupath/new_detection_filter_invariants.groovy'
QUPATH = Path('/home/smith6jt/QuPath/bin/QuPath')
REQUIRED_STEPS = tuple(s.id for stage in AU.STAGES for s in AU.stage_steps(stage))
_HASH_CACHE: dict[tuple, str] = {}


def read_json(path: Path):
    return json.loads(Path(path).read_text())


def digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def file_stat(path: Path) -> dict:
    p = Path(path).resolve(); s = p.stat()
    return {'path': str(p), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns,
            'ctime_ns': s.st_ctime_ns, 'inode': s.st_ino}


def sha(path: Path) -> str:
    before = file_stat(path); key = tuple(before.values())
    if key not in _HASH_CACHE:
        h = hashlib.sha256()
        with Path(path).open('rb') as stream:
            for chunk in iter(lambda: stream.read(1 << 20), b''): h.update(chunk)
        if file_stat(path) != before: raise ValueError(f'File changed while hashing: {path}')
        _HASH_CACHE[key] = h.hexdigest()
    return _HASH_CACHE[key]


def signature(path: Path) -> dict:
    before = file_stat(path); result = {**before, 'sha256': sha(path)}
    if file_stat(path) != before: raise ValueError(f'File changed while signing: {path}')
    return result


def record(path: Path) -> dict:
    s = signature(path)
    return {k: s[k] for k in ('path', 'size', 'sha256')}


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, prefix=path.name + '.', delete=False) as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
        temporary = Path(f.name)
    os.replace(temporary, path)


LOADED_CODE = {str(p): sha(p) for p in (Path(__file__), SCRIPT, HELPER)}


def code_fingerprint() -> dict:
    current = {str(p): sha(Path(p)) for p in LOADED_CODE}
    if current != LOADED_CODE: raise ValueError('Filter code changed since import; restart before processing')
    return current


def unique_ids(values, label: str) -> set[str]:
    seq = list(values)
    if any(not isinstance(v, str) or not v.strip() for v in seq) or len(set(seq)) != len(seq):
        raise ValueError(f'{label}: null, empty, non-string or duplicate ID')
    return set(seq)


def source_seal(source: Path, tag: str) -> dict:
    """Seal a fully completed raw image, including every checkpoint output and raw evidence file."""
    source = Path(source).resolve()
    area = source / 'augmentation' / tag
    state_path = area / 'state.json'; state = read_json(state_path)
    if (state.get('schema') != 'uniform.tissue_augmentation/1' or state.get('in_progress')
            or set(REQUIRED_STEPS) - set(state.get('steps', {}))):
        raise ValueError(f'{tag}: raw augmentation is incomplete')
    state_hash = sha(state_path)
    checkpoint_outputs = {}
    for step in REQUIRED_STEPS:
        item = state['steps'][step]
        if not item.get('fingerprint') or not isinstance(item.get('outputs'), dict):
            raise ValueError(f'{tag}: malformed raw checkpoint {step}')
        for name, expected in item['outputs'].items():
            p = Path(name)
            if sha(p) != expected: raise ValueError(f'{tag}: stale raw checkpoint output {p}')
            checkpoint_outputs[str(p.resolve())] = expected
    entry = CT.manifest(source)[tag]
    paths = {
        'project_manifest': source / 'project_manifest.json', 'project': source / 'qupath/project.qpproj',
        'server': source / 'qupath/data' / str(entry['entry_id']) / 'server.json', 'qpdata': AU.qpdata(source, tag),
        'cells': source / 'cells' / f'{tag}.parquet', 'final_cells': source / 'final_cells' / f'{tag}.parquet',
        'context': source / 'context' / f'{tag}.parquet', 'geometry': source / 'cell_geojson' / f'{tag}.geojson',
        'annotations': source / 'annotations' / f'{tag}.tsv', 'islets': source / 'islets' / f'{tag}.geojson',
        'tissue': source / 'inputs/tissue' / f'{tag}.geojson', 'augmentation_state': state_path,
        'verification': area / 'verification.json', 'original_inventory': area / 'original_inventory.json',
        'raw_new_ids': area / 'new_object_ids.json',
    }
    files = {k: record(p) for k, p in paths.items()}
    for key in ('cells', 'final_cells', 'context', 'geometry', 'annotations', 'islets', 'tissue',
                'verification', 'original_inventory', 'raw_new_ids'):
        item = files[key]
        if checkpoint_outputs.get(item['path']) != item['sha256']:
            raise ValueError(f'{tag}: {key} is not pinned by a completed raw checkpoint')
    if files['qpdata']['sha256'] != state['last_qpdata_sha256']:
        raise ValueError(f'{tag}: source qpdata changed after raw augmentation')
    if files['tissue']['sha256'] != state['inputs']['outline_sha256']:
        raise ValueError(f'{tag}: source Tissue differs from accepted augmentation outline')
    original = unique_ids(read_json(paths['original_inventory'])['cells'], 'original IDs')
    new = unique_ids(read_json(paths['raw_new_ids']), 'raw new IDs')
    report = read_json(paths['verification'])
    if original & new or any(report.get(k) != n for k, n in {
        'original_cells': len(original), 'new_cells': len(new), 'total_cells': len(original | new)}.items()):
        raise ValueError(f'{tag}: invalid raw inventory partition/counts')
    evidence_paths = set(area.glob('*'))
    evidence_paths.update((source / 'provenance' / tag).glob('augment_*.json'))
    evidence_paths.update((source / 'logs' / tag).glob('augment_*.log'))
    evidence = {str(p.resolve()): record(p) for p in sorted(evidence_paths) if p.is_file()}
    result = {'schema': 'uniform.raw_augmentation_seal/1', 'source_run': str(source), 'image': tag,
              'entry': entry, 'channels': CT.channels(source, entry['entry_id']),
              'source_image': file_stat(Path(entry['image_path'])), 'files': files,
              'checkpoint_outputs': checkpoint_outputs, 'preserved_source_evidence': evidence,
              'counts': {'original_cells': len(original), 'raw_new_cells': len(new), 'total_cells': len(original | new)},
              'original_ids_sha256': digest(sorted(original)), 'raw_new_ids_sha256': digest(sorted(new))}
    if sha(state_path) != state_hash or result['files']['augmentation_state']['sha256'] != state_hash:
        raise ValueError(f'{tag}: raw augmentation changed during seal')
    return result


def decision_projection(document: dict, tag: str) -> dict:
    if document.get('schema') != DECISIONS_SCHEMA: raise ValueError('Unsupported decisions schema')
    rule = document.get('rule', {})
    if (any(not isinstance(rule.get(k), str) or not rule[k].strip() for k in ('id', 'version', 'decision_basis'))
            or not isinstance(rule.get('parameters'), dict)):
        raise ValueError('A common explicit rule id/version/parameters/decision_basis is required')
    return {'schema': DECISIONS_SCHEMA, 'rule': rule, 'image': tag, 'record': document['images'][tag]}


def validate_decisions(document: dict, tag: str, seal: dict) -> tuple[dict, list[str], list[str], list[str]]:
    projected = decision_projection(document, tag); item = projected['record']
    if item.get('source_seal') != seal: raise ValueError(f'{tag}: decisions source seal is stale')
    review = item.get('review', {})
    if review.get('status') != 'reviewed' or not review.get('reviewer') or not review.get('evidence'):
        raise ValueError(f'{tag}: reviewed scientific evidence is required')
    for evidence in review['evidence']:
        if sha(Path(evidence['path'])) != evidence.get('sha256'):
            raise ValueError(f'{tag}: reviewed evidence changed')
    raw = unique_ids(read_json(Path(seal['files']['raw_new_ids']['path'])), 'raw new IDs')
    decisions = item.get('decisions', [])
    ids = unique_ids([d['object_id'] for d in decisions], 'decision IDs')
    if ids != raw: raise ValueError(f'{tag}: decisions must cover exactly every raw new detection')
    accepted, rejected, deferred = [], [], []
    for d in decisions:
        if d.get('decision') not in ('accept', 'reject', 'defer') or not isinstance(d.get('reason'), str) or not d['reason'].strip():
            raise ValueError(f'{tag}: every new ID requires explicit accept/reject/defer and reason')
        {'accept':accepted, 'reject':rejected, 'defer':deferred}[d['decision']].append(d['object_id'])
    return projected, sorted(accepted), sorted(rejected), sorted(deferred)


def _features(path: Path) -> tuple[dict | list, list]:
    value = read_json(path)
    return value, value['features'] if isinstance(value, dict) else value


def filter_artifacts(source: Path, stage: Path, tag: str, retained: set[str], rejected: set[str]) -> None:
    """Subset complete rows/features; never recalculate or coerce retained quantitative values."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    for folder in ('cells', 'final_cells', 'context'):
        src = source / folder / f'{tag}.parquet'; dst = stage / folder / src.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        table = pq.read_table(src); ids = table['Object ID'].to_pylist()
        actual = unique_ids(ids, folder)
        if actual != retained | rejected: raise ValueError(f'{folder}: source inventory mismatch')
        if not rejected: shutil.copy2(src, dst)
        else: pq.write_table(table.filter(pa.array([x in retained for x in ids])), dst)
        evidence = stage / 'filter' / tag / 'excluded' / f'{folder}.parquet'; evidence.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table.filter(pa.array([x in rejected for x in ids])), evidence)
    src = source / 'cell_geojson' / f'{tag}.geojson'; dst = stage / 'cell_geojson' / src.name
    value, features = _features(src); ids = unique_ids([f['id'] for f in features], 'GeoJSON')
    if ids != retained | rejected: raise ValueError('Source GeoJSON inventory differs')
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not rejected: shutil.copy2(src, dst)
    else:
        filtered = [f for f in features if f['id'] in retained]
        atomic_json(dst, {**value, 'features': filtered} if isinstance(value, dict) else filtered)
    atomic_json(stage / 'filter' / tag / 'excluded/geometry.geojson',
                {'type': 'FeatureCollection', 'features': [f for f in features if f['id'] in rejected]})


def verify_artifacts(source: Path, stage: Path, tag: str, retained: set[str], rejected: set[str]) -> dict:
    import numpy as np
    import pyarrow.parquet as pq
    for folder in ('cells', 'final_cells', 'context'):
        old = pq.read_table(source / folder / f'{tag}.parquet'); new = pq.read_table(stage / folder / f'{tag}.parquet')
        if not old.schema.equals(new.schema, check_metadata=True): raise ValueError(f'{folder}: schema changed')
        oi = old['Object ID'].to_pylist(); ni = new['Object ID'].to_pylist()
        if unique_ids(ni, folder) != retained or ni != [x for x in oi if x in retained]:
            raise ValueError(f'{folder}: filtered inventory/order differs')
        keep = np.array([x in retained for x in oi])
        for name in old.column_names:
            a = old[name].combine_chunks(); b = new[name].combine_chunks()
            if a.type != b.type: raise ValueError(f'{folder}/{name}: dtype changed')
            if not a.is_valid().filter(__import__('pyarrow').array(keep)).equals(b.is_valid()):
                raise ValueError(f'{folder}/{name}: retained null mask changed')
            # Compare Arrow scalar values, preserving float payloads (including NaN and signed zero).
            left = a.to_pandas().to_numpy()[keep]; right = b.to_pandas().to_numpy()
            if left.dtype.kind in 'fc':
                equal = left.dtype == right.dtype and left.tobytes() == right.tobytes()
            else:
                equal = old.select([name]).filter(__import__('pyarrow').array(keep)).equals(new.select([name]))
            if not equal: raise ValueError(f'{folder}/{name}: retained value changed')
    _, original = _features(source / 'cell_geojson' / f'{tag}.geojson')
    _, current = _features(stage / 'cell_geojson' / f'{tag}.geojson')
    if current != [f for f in original if f['id'] in retained]: raise ValueError('Retained cell/nucleus GeoJSON changed')
    for folder, suffix in (('islets', 'geojson'), ('inputs/tissue', 'geojson')):
        if sha(source / folder / f'{tag}.{suffix}') != sha(stage / folder / f'{tag}.{suffix}'):
            raise ValueError(f'{folder}: annotation/Tissue bytes changed')
    if not rejected:
        for folder, suffix in (('cells','parquet'), ('final_cells','parquet'), ('context','parquet'),
                               ('cell_geojson','geojson'), ('annotations','tsv')):
            if sha(source/folder/f'{tag}.{suffix}') != sha(stage/folder/f'{tag}.{suffix}'):
                raise ValueError(f'{folder}: no-op generation must be byte-identical')
        if sha(AU.qpdata(source, tag)) != sha(AU.qpdata(stage, tag)):
            raise ValueError('No-op qpdata is not byte-identical')
    return {'retained_cells': len(retained), 'excluded_cells': len(rejected),
            'retained_rows_and_geometry_exact': True, 'endocrine_and_tissue_bytes_unchanged': True,
            'noop_byte_identical': not rejected}


def verify_annotations(source: Path, stage: Path, tag: str, rejected: set[str], native: dict) -> None:
    import pandas as pd
    src = pd.read_csv(source/'annotations'/f'{tag}.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('Object ID')
    dst = pd.read_csv(stage/'annotations'/f'{tag}.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('Object ID')
    if not src.index.is_unique or not dst.index.is_unique or set(src.index) != set(dst.index) or set(src.columns) != set(dst.columns):
        raise ValueError('Annotation table inventory/schema changed')
    columns = [c for c in src if c != 'Num Detections']
    if not src[columns].equals(dst.loc[src.index, columns]):
        raise ValueError('Stored annotation measurements/metadata changed')
    if 'Num Detections' in src:
        old = pd.to_numeric(src['Num Detections']); new = pd.to_numeric(dst.loc[src.index, 'Num Detections'])
        if ((old-new) < 0).any() or ((old-new) > len(rejected)).any():
            raise ValueError('Dynamic annotation detection counts changed inconsistently')
        # Parent ancestry is separately fingerprinted. QuPath's exported dynamic
        # count follows its measurement table's spatial-membership semantics.
        expected=native['annotation_export_detection_counts']
        if set(expected) != set(new.index):
            raise ValueError('Native annotation export-count inventory differs')
        if any(float(new.loc[oid]) != expected[oid] for oid in new.index):
            raise ValueError('Annotation export dynamic counts differ from reloaded native hierarchy')


def run_native(stage: Path, tag: str, original: set[str], raw_new: set[str], rejected: set[str]) -> dict:
    """Two fresh QuPath processes: deletion/save, then read-only reload verification."""
    area = stage/'filter'/tag; config = area/'native_config.json'
    common = {'helper': str(HELPER), 'original_ids': sorted(original), 'raw_new_ids': sorted(raw_new),
              'excluded_ids': sorted(rejected), 'before_inventory': str(area/'native_before.json'),
              'annotation_output': str(stage/'annotations'/f'{tag}.tsv')}
    for mode in ('apply', 'verify'):
        report = area/f'native_{mode}.json'
        atomic_json(config, {**common, 'mode': mode, 'report': str(report)})
        log = area/f'native_{mode}.log'
        with log.open('w') as stream:
            result = subprocess.run([str(QUPATH), 'script', '-p', str(stage/'qupath/project.qpproj'), '-i', tag,
                                     '-a', str(config), str(SCRIPT)], stdin=subprocess.DEVNULL, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=1800,
                                    env={**os.environ, 'JAVA_TOOL_OPTIONS':'-Xmx64g -XX:ActiveProcessorCount=8'})
        if result.returncode or not report.is_file(): raise RuntimeError(f'{tag}: native {mode} failed; {log}')
        value = read_json(report)
        if (not value.get('ok') or not value.get('retained_full_signatures_match')
                or not value.get('annotation_full_signatures_match') or value.get('excluded_cells') != len(rejected)):
            raise ValueError(f'{tag}: incomplete native {mode} verification')
    return read_json(area/'native_verify.json')


def _copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(src, dst)
    AU.check_independent(src, dst)


def _output_records(stage: Path, tag: str, entry_id: int) -> dict:
    keys = {'cells':f'cells/{tag}.parquet','final_cells':f'final_cells/{tag}.parquet','context':f'context/{tag}.parquet',
            'geometry':f'cell_geojson/{tag}.geojson','annotations':f'annotations/{tag}.tsv','islets':f'islets/{tag}.geojson',
            'tissue':f'inputs/tissue/{tag}.geojson','qpdata':f'qupath/data/{entry_id}/data.qpdata',
            'server':f'qupath/data/{entry_id}/server.json','project':'qupath/project.qpproj','project_manifest':'project_manifest.json'}
    for p in (stage/'filter'/tag).rglob('*'):
        if p.is_file() and p.name != 'manifest.json': keys['evidence:'+str(p.relative_to(stage))] = str(p.relative_to(stage))
    return {key:{'relative_path':name,'sha256':sha(stage/name),'size':(stage/name).stat().st_size} for key,name in keys.items()}


def process(source: Path, out: Path, decisions_path: Path, tag: str) -> dict:
    source, out, decisions_path = map(lambda p: Path(p).resolve(), (source,out,decisions_path))
    if out == source or source in out.parents or out in source.parents: raise ValueError('Source and QC roots must be independent siblings')
    seal = source_seal(source,tag); projection, accepted, rejected, deferred = validate_decisions(read_json(decisions_path),tag,seal)
    expected_code = code_fingerprint(); completion = out/'filter'/tag/'manifest.json'
    if completion.exists():
        current = accepted_inputs(out,tag)
        if current['decision_projection_sha256'] != digest(projection): raise ValueError('Completed decisions changed; use a new generation')
        return read_json(completion)
    out.mkdir(parents=True,exist_ok=True)
    lock_path=out/'filter'/f'{tag}.lock';lock_path.parent.mkdir(parents=True,exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        if completion.exists():
            accepted_inputs(out,tag)
            return read_json(completion)
        intent={'source_seal_sha256':digest(seal),'decision_projection_sha256':digest(projection),
                'common_rule_sha256':digest(projection['rule']),'code_sha256':expected_code}
        pending=out/'filter'/tag/'pending.json'
        generation=out/'filter/generation.json'
        if generation.exists() and read_json(generation)!={'schema':SCHEMA,'source_run':str(source),'common_rule_sha256':intent['common_rule_sha256']}:
            raise ValueError('QC generation source/common rule differs')
        if pending.exists():
            if read_json(pending)!=intent:raise ValueError('Pending filter attempt has different inputs; use a new generation')
        else:
            for name in ('cells','final_cells','context','geometry','annotations','islets','tissue','qpdata','server'):
                target=out/Path(seal['files'][name]['path']).relative_to(source)
                if target.exists():raise ValueError(f'Refusing to overwrite an unowned existing output: {target}')
            atomic_json(pending,intent)
        original=unique_ids(read_json(Path(seal['files']['original_inventory']['path']))['cells'],'original IDs')
        excluded=set(rejected)|set(deferred);raw_new=set(accepted)|excluded;retained=original|set(accepted)
        staging=out/'filter'/'.staging';staging.mkdir(exist_ok=True)
        # Each attempt is preserved if it fails; retries start from the sealed raw source, never a partial save.
        stage=Path(tempfile.mkdtemp(prefix=tag+'.',dir=staging));area=stage/'filter'/tag;area.mkdir(parents=True)
        atomic_json(area/'source_seal.json',seal);atomic_json(area/'decisions.json',projection)
        for name,ids in [('raw_new_ids',raw_new),('accepted_new_ids',accepted),('rejected_new_ids',rejected),('deferred_new_ids',deferred)]:
            atomic_json(area/(name+'.json'),sorted(ids))
        for key in ('project_manifest','project','server','qpdata','annotations','islets','tissue'):
            src=Path(seal['files'][key]['path']);_copy(src,stage/src.relative_to(source))
        filter_artifacts(source,stage,tag,retained,excluded)
        native=run_native(stage,tag,original,raw_new,excluded)
        checks=verify_artifacts(source,stage,tag,retained,excluded);verify_annotations(source,stage,tag,excluded,native)
        if source_seal(source,tag)!=seal or code_fingerprint()!=expected_code:
            raise ValueError('Source/code changed during filtering; output not published')
        if digest(decision_projection(read_json(decisions_path),tag))!=digest(projection):
            raise ValueError('Decisions changed during filtering; output not published')
        atomic_json(area/'verification.json',{**checks,'native':native})
        outputs=_output_records(stage,tag,seal['entry']['entry_id'])
        result={'schema':SCHEMA,'image':tag,'output_run':str(out),'source_run':str(source),
                'source_seal':seal,'decision_source':str(decisions_path),'decision_projection_sha256':digest(projection),
                'common_rule_sha256':digest(projection['rule']),'code_sha256':expected_code,
                'original_ids_sha256':seal['original_ids_sha256'],'raw_new_ids':sorted(raw_new),
                'accepted_new_ids':accepted,'rejected_new_ids':rejected,'deferred_new_ids':deferred,
                'counts':{'original_cells':len(original),'raw_new_cells':len(raw_new),'accepted_new_cells':len(accepted),
                          'rejected_new_cells':len(rejected),'deferred_new_cells':len(deferred),'excluded_new_cells':len(excluded),'total_cells':len(retained)},
                'outputs':outputs,'checks':checks,'status':'complete_requires_current_geometry_review',
                'completed_at_utc':datetime.now(timezone.utc).isoformat()}
        # Shared metadata must remain byte-identical. A root lock serializes short publication only.
        with (out/'filter/publish.lock').open('a') as publication:
            fcntl.flock(publication,fcntl.LOCK_EX)
            generation=out/'filter/generation.json'
            root_record={'schema':SCHEMA,'source_run':str(source),'common_rule_sha256':result['common_rule_sha256']}
            if generation.exists() and read_json(generation)!=root_record: raise ValueError('QC generation source/common rule differs')
            atomic_json(generation,root_record)
            for item in outputs.values():
                rel=item['relative_path'];src=stage/rel;dest=out/rel;dest.parent.mkdir(parents=True,exist_ok=True)
                if dest.exists() and rel in ('project_manifest.json','qupath/project.qpproj'):
                    if sha(dest)!=item['sha256']: raise ValueError('Shared project metadata changed')
                    continue
                temporary=dest.with_name(dest.name+'.filter-part');_copy(src,temporary);os.replace(temporary,dest)
            # Readers never trust published files without this final atomic marker.
            if source_seal(source,tag)!=seal: raise ValueError('Source changed while publishing; completion withheld')
            atomic_json(completion,result)
        accepted_inputs(out,tag)
        return result


def accepted_inputs(out: Path, tag: str) -> dict | None:
    """Validated current filtered input snapshot; scientific review/release is a separate gate."""
    out=Path(out).resolve();path=out/'filter'/tag/'manifest.json'
    if not path.is_file():return None
    initial=sha(path);m=read_json(path)
    if m.get('schema')!=SCHEMA or m.get('output_run')!=str(out) or m.get('image')!=tag:
        raise ValueError('Unsupported or misplaced filter manifest')
    if m['code_sha256']!=code_fingerprint():raise ValueError('Filter implementation changed since completion')
    seal=source_seal(Path(m['source_run']),tag)
    if seal!=m['source_seal']:raise ValueError(f'{tag}: immutable raw source changed')
    projection,accepted,rejected,deferred=validate_decisions(read_json(Path(m['decision_source'])),tag,seal)
    if digest(projection)!=m['decision_projection_sha256'] or accepted!=m['accepted_new_ids'] or rejected!=m['rejected_new_ids'] or deferred!=m['deferred_new_ids']:
        raise ValueError(f'{tag}: reviewed decisions changed')
    if digest(projection['rule'])!=m['common_rule_sha256']:raise ValueError('Common rule changed')
    expected_counts={'original_cells':seal['counts']['original_cells'],'raw_new_cells':seal['counts']['raw_new_cells'],
                     'accepted_new_cells':len(accepted),'rejected_new_cells':len(rejected),'deferred_new_cells':len(deferred),
                     'excluded_new_cells':len(rejected)+len(deferred),'total_cells':seal['counts']['original_cells']+len(accepted)}
    if m['counts']!=expected_counts or m['raw_new_ids']!=sorted(accepted+rejected+deferred):
        raise ValueError('Filtered inventory partition/counts changed')
    if m['status']!='complete_requires_current_geometry_review' or not m['checks'].get('retained_rows_and_geometry_exact'):
        raise ValueError('Filter completion is not verified')
    generation=read_json(out/'filter/generation.json')
    if generation!={'schema':SCHEMA,'source_run':m['source_run'],'common_rule_sha256':m['common_rule_sha256']}:
        raise ValueError('Generation lineage changed')
    for key,item in m['outputs'].items():
        p=(out/item['relative_path']).resolve()
        if not p.is_relative_to(out) or sha(p)!=item['sha256'] or p.stat().st_size!=item['size']:
            raise ValueError(f'{tag}: filtered output changed: {key}')
    if sha(path)!=initial:raise ValueError('Filter manifest changed during validation')
    files={k:signature(out/m['outputs'][k]['relative_path']) for k in
           ('cells','final_cells','context','geometry','annotations','islets','tissue','qpdata','server')}
    files['verification']=signature(out/'filter'/tag/'verification.json')
    return {'run':str(out),'image':tag,'generation':SCHEMA,'files':files,'entry':seal['entry'],
            'channels':seal['channels'],'source_image':seal['source_image'],'total_cells':m['counts']['total_cells'],
            'counts':m['counts'],'filter_manifest_sha256':initial,'filter_manifest_path':str(path),
            'decision_projection_sha256':m['decision_projection_sha256'],
            'source_augmentation_state_sha256':seal['files']['augmentation_state']['sha256'],
            'source_run':m['source_run']}


def _structure_inputs(source: Path, tag: str, seal: dict) -> dict:
    import redraw_structures as R
    row=next(r for r in CT.cohort_rows() if r['image_name']==tag)
    return {'tissue_sha256':seal['files']['tissue']['sha256'],'qpdata_sha256':seal['files']['qpdata']['sha256'],
            'script_sha256':R.sha(R.SCRIPT),'driver_sha256':R.sha(Path(R.__file__)),
            'threshold_code_sha256':R.sha(Path(R.ST.__file__)),'threshold_inputs':R.threshold_inputs(source,row),
            'classifiers':{c:R.sha(R.CLASSIFIERS/f'{c}.json') for c in ('Nerve','Capillary','Lymphatics')}}


def certify_structures(out: Path, tag: str) -> dict:
    """Copy established masks with explicit unchanged-input lineage; never retag an old redraw record."""
    out=Path(out).resolve();snapshot=accepted_inputs(out,tag)
    if snapshot is None:raise ValueError('Filtered generation is not complete')
    existing=out/'filter'/tag/'structures_reuse.json'
    if existing.is_file():
        structure_reuse_current(out,tag)
        return read_json(existing)
    source=Path(snapshot['source_run']);m=read_json(out/'filter'/tag/'manifest.json');seal=m['source_seal']
    old_path=source/'structures/redraw_records'/f'{tag}.json';old_hash=sha(old_path);old=read_json(old_path)
    expected=_structure_inputs(source,tag,seal)
    if old['inputs']!=expected or set(old['classes'])!={'Nerve','Capillary','Lymphatic'}:
        raise ValueError('Source structure inference is stale/incomplete')
    outputs={}
    for cls,item in old['classes'].items():
        src=source/'structures'/cls/f'{cls}__{tag}.geojson';dst=out/'structures'/cls/src.name
        if sha(src)!=item['sha256'] or item.get('edge_touching_retained') is not True:raise ValueError('Source structure mask changed')
        if dst.exists() and sha(dst)!=item['sha256']:raise ValueError('Existing QC structure differs; refusing overwrite')
        if not dst.exists():_copy(src,dst)
        outputs[cls]={'relative_path':str(dst.relative_to(out)),'source':record(src),'sha256':sha(dst)}
    if accepted_inputs(out,tag)!=snapshot or sha(old_path)!=old_hash or _structure_inputs(source,tag,seal)!=expected:
        raise ValueError('Inputs changed during structure reuse')
    cert={'schema':'uniform.established_structure_filter_reuse/1','image':tag,'filter_manifest_sha256':snapshot['filter_manifest_sha256'],
          'source_record':record(old_path),'source_inference_inputs':expected,'outputs':outputs,
          'current_qpdata_sha256':snapshot['files']['qpdata']['sha256'],
          'current_tissue_sha256':snapshot['files']['tissue']['sha256'],
          'proof':'Fresh native reload verified all annotation ROIs/metadata and all retained cells unchanged; only rejected leaf detections removed. Pixel inference inputs are unchanged.',
          'created_at_utc':datetime.now(timezone.utc).isoformat()}
    atomic_json(out/'filter'/tag/'structures_reuse.json',cert)
    structure_reuse_current(out,tag)
    return cert


def structure_reuse_current(out: Path, tag: str) -> dict | None:
    out=Path(out).resolve();path=out/'filter'/tag/'structures_reuse.json'
    if not path.is_file():return None
    snap=accepted_inputs(out,tag);c=read_json(path)
    if snap is None or c['filter_manifest_sha256']!=snap['filter_manifest_sha256']:
        raise ValueError('Structure reuse generation changed')
    if c['current_qpdata_sha256']!=snap['files']['qpdata']['sha256'] or c['current_tissue_sha256']!=snap['files']['tissue']['sha256']:
        raise ValueError('Structure reuse geometry changed')
    if sha(Path(c['source_record']['path']))!=c['source_record']['sha256']:raise ValueError('Source structure record changed')
    manifest=read_json(out/'filter'/tag/'manifest.json')
    if _structure_inputs(Path(snap['source_run']),tag,manifest['source_seal'])!=c['source_inference_inputs']:
        raise ValueError('Source structure classifier/threshold/code inputs changed')
    for item in c['outputs'].values():
        if sha(Path(item['source']['path']))!=item['sha256'] or sha(out/item['relative_path'])!=item['sha256']:
            raise ValueError('Reused structure mask changed')
    return {'certificate':str(path),'certificate_sha256':sha(path),'outputs':c['outputs']}


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run',type=Path);p.add_argument('--out',type=Path);p.add_argument('--decisions',type=Path)
    p.add_argument('--images',nargs='+',required=True);p.add_argument('--dry-run',action='store_true')
    p.add_argument('--seal-only',action='store_true');p.add_argument('--verify-only',action='store_true')
    p.add_argument('--certify-structures',action='store_true')
    a=p.parse_args()
    if a.seal_only:
        if a.source_run is None:p.error('--seal-only requires --source-run')
        print(json.dumps({t:source_seal(a.source_run,t) for t in a.images},indent=2));return 0
    if a.out is None:p.error('--out is required')
    if a.verify_only or a.certify_structures:
        for tag in a.images:
            result=certify_structures(a.out,tag) if a.certify_structures else accepted_inputs(a.out,tag)
            if result is None:raise ValueError(tag+': filtered generation is incomplete')
            print(tag+': '+('structures certified' if a.certify_structures else 'filtered inputs current'))
        return 0
    if a.source_run is None or a.decisions is None:p.error('--source-run and --decisions are required')
    for tag in a.images:
        if a.dry_run:
            seal=source_seal(a.source_run,tag);_,accepted,rejected,deferred=validate_decisions(read_json(a.decisions),tag,seal)
            print(json.dumps({'image':tag,'accepted_new':len(accepted),'rejected_new':len(rejected),'deferred_new':len(deferred),'mutation':False}))
        else:
            result=process(a.source_run,a.out,a.decisions,tag)
            print(json.dumps({'image':tag,**result['counts'],'status':result['status']}))
    return 0


if __name__=='__main__':raise SystemExit(main())
