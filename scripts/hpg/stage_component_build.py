#!/usr/bin/env python3
"""Prepare fresh, hash-pinned inputs for a diagnostic endocrine component build.

Capture runs on the workstation; materialize and prepare run through the sealed
HPG site wrapper. No command submits a job, extracts pixels, assigns identities,
accepts biology, creates release attestations or starts a service.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
COPIED_DIRS = ('typing/context', 'typing/identity', 'typing/coverage', 'cell_structure_distances')
COPIED_FILES = ('structure_generation.json', 'structure_metrics_manifest.json', 'structure_region_metrics.csv')
SHARED_DIRS = ('cells', 'cell_geojson', 'context', 'qupath', 'inputs', 'islets', 'structures', 'annotations',
               'typing/redsea')
MANIFEST = 'initial_source_package.json'


def sha(path):
    path = Path(path)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b''):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino) != (
            after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino):
        raise ValueError(f'File changed while hashing: {path}')
    return digest.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def independent_output(output, sources):
    """Reject source nesting before any mkdir, copy, or link operation."""
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    for source in sources:
        if output == source or output.is_relative_to(source) or source.is_relative_to(output):
            raise ValueError(f'Output must be independent of protected source {source}: {output}')


def inventory(root, folders, files=()):
    result = set(files)
    for folder in folders:
        source = root / folder
        if not source.is_dir():
            raise FileNotFoundError(source)
        result.update(str(p.relative_to(root)) for p in source.rglob('*') if p.is_file())
    return sorted(result)


def checked(root, relative, record):
    if Path(relative).is_absolute() or '..' in Path(relative).parts:
        raise ValueError(f'Unsafe source relative path: {relative}')
    path = root / relative
    if path.stat().st_size != record['bytes'] or sha(path) != record['sha256']:
        raise ValueError(f'Source hash/size mismatch: {path}')
    return path


def copy_checked(source, target, digest):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    shutil.copy2(source, target)
    if sha(target) != digest or sha(source) != digest:
        raise ValueError(f'Source changed while copying: {source}')


def capture(args):
    source, output = args.source.resolve(strict=True), args.output.resolve()
    independent_output(output, (source,))
    output.mkdir(parents=True, exist_ok=False)
    copied = inventory(source, COPIED_DIRS, COPIED_FILES)
    copied += [str(p.relative_to(source)) for p in sorted((source / 'typing').glob('*.csv'))]
    shared = inventory(source, SHARED_DIRS, ('project_manifest.json',))
    tasks = [(relative, 'copied') for relative in copied] + [(relative, 'shared') for relative in shared]
    def one(task):
        relative, kind = task
        path = source / relative
        record = {'bytes': path.stat().st_size, 'sha256': sha(path)}
        if kind == 'copied':
            copy_checked(path, output / relative, record['sha256'])
            (output / relative).chmod(0o444)
        return relative, kind, record
    records = {'copied': {}, 'shared': {}}
    with ThreadPoolExecutor(args.hash_workers) as pool:
        for relative, kind, record in pool.map(one, tasks):
            records[kind][relative] = record
    record = {'schema': 'uniform.initial_component_source_package/1',
              'created_utc': datetime.now(timezone.utc).isoformat(), 'source_run': str(source),
              'scientific_baseline': args.baseline_description,
              'historical_generation_preserved_without_rebinding': True, **records}
    write_json(output / MANIFEST, record)
    (output / MANIFEST).chmod(0o444)
    print(json.dumps({'package': str(output), 'manifest_sha256': sha(output / MANIFEST),
                      'copied_files': len(copied), 'copied_bytes': sum(r['bytes'] for r in records['copied'].values()),
                      'shared_files': len(shared), 'shared_bytes': sum(r['bytes'] for r in records['shared'].values())}), flush=True)


def package_record(package):
    record = json.loads((package / MANIFEST).read_text())
    if record.get('schema') != 'uniform.initial_component_source_package/1':
        raise ValueError('Unsupported initial source package')
    if set(record['copied']) & set(record['shared']):
        raise ValueError('Overlapping copied and shared package inputs')
    return record


def materialize(args):
    package, shared, output = args.package.resolve(strict=True), args.shared_run.resolve(strict=True), args.output.resolve()
    record = package_record(package)
    independent_output(output, (package, shared))
    # Validate every shared byte before creating any live links.
    tasks = [(package, p, r) for p, r in record['copied'].items()] + [
             (shared, p, r) for p, r in record['shared'].items()]
    with ThreadPoolExecutor(args.hash_workers) as pool:
        list(pool.map(lambda t: checked(*t), tasks))
    output.mkdir(parents=True, exist_ok=False)
    for relative, item in record['copied'].items():
        copy_checked(package / relative, output / relative, item['sha256'])
        (output / relative).chmod(0o444)
    for relative in record['shared']:
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to((shared / relative).resolve(strict=True))
    write_json(output / 'materialized_source.json', {
        'schema': 'uniform.materialized_initial_component_source/1', 'run': str(output),
        'source_package': str(package), 'package_manifest_sha256': sha(package / MANIFEST),
        'shared_run': str(shared), 'source_manifest': record,
        'note': 'Historical generation bytes retained. This is a calibration/identity source, not a newly validated geometry run.'})
    print(f'Materialized exact initial source {output}', flush=True)


def prepare(args):
    sys.path.insert(0, str(ROOT / 'scripts/uniform/typing'))
    sys.path.insert(0, str(ROOT / 'scripts/uniform'))
    import pandas as pd
    import area_coverage as AC
    import calibration_context as CC
    import component_selection as CS
    import codex_typing as CT
    import structure_generation as SG
    import structure_metrics as SM
    initial, endocrine, structures = (p.resolve(strict=True) for p in (
        args.nonendocrine_source, args.endocrine_source, args.structure_source))
    run = args.run.resolve()
    independent_output(run, (initial, endocrine, structures))
    materialized = json.loads((initial / 'materialized_source.json').read_text())
    if materialized.get('schema') != 'uniform.materialized_initial_component_source/1' or materialized.get('run') != str(initial):
        raise ValueError('Initial source must be the explicitly materialized, verified package')
    package = Path(materialized['source_package'])
    if sha(package / MANIFEST) != materialized['package_manifest_sha256']:
        raise ValueError('Initial source package manifest changed')
    source = materialized['source_manifest']
    for group in ('copied', 'shared'):
        for relative, item in source[group].items():
            checked(initial, relative, item)
    if args.classification_control and structures != initial:
        raise ValueError('Classification control must explicitly use the initial source structures/context')
    if not args.classification_control:
        SG.validate(structures)
    run.mkdir(parents=True, exist_ok=False)
    # Shared files are strictly inputs. All output folders below are separate copies.
    for relative, item in source['shared'].items():
        path = checked(structures, relative, item)
        target = run / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(path.resolve(strict=True))
    for relative in inventory(structures, ('typing/context', 'cell_structure_distances'),
                              ('structure_metrics_manifest.json', 'structure_region_metrics.csv')):
        copy_checked(structures / relative, run / relative, sha(structures / relative))
    generation = structures / 'structure_generation.json'
    provenance = run / 'source_provenance/structure_generation.json'
    copy_checked(generation, provenance, sha(generation))
    if args.classification_control:
        copy_checked(generation, run / 'structure_generation.json', sha(generation))
    else:
        SG.pin_derivation(run, structures, selection_reason=args.selection_reason)
        # Includes the unchanged protected Collagen IV validator; no provenance
        # relabeling or waived checks are needed when immutable links resolve to source.
        SM.validate_manifest(run, [row['image_name'] for row in CT.cohort_rows()])
    old, other = CS.typing_root(endocrine), initial / 'typing'
    for name in CS.MIXED:
        a, b = (pd.read_csv(p / name, float_precision='round_trip') for p in (old, other))
        CS.mixed_table(a, b).to_csv(run / 'typing' / name, index=False)
    for name in CS.UNCHANGED:
        copy_checked(other / name, run / 'typing' / name, sha(other / name))
    tags = [row['image_name'] for row in CT.cohort_rows()]
    if len(tags) != 20 or len(set(tags)) != 20:
        raise ValueError('Expected the explicit 20-image cohort')
    if any(str(row['case_id']).lstrip('0') == '78' for row in CT.cohort_rows()):
        raise ValueError('Excluded donor 78 cannot enter a component build')
    derivations = {}
    for tag in tags:
        frames, origins = [], []
        for root in (old, other):
            path = root / 'coverage' / f'{tag}.parquet'
            sidecar = path.with_suffix('.json')
            metadata = json.loads(sidecar.read_text())
            if metadata['parquet_sha256'] != sha(path):
                raise ValueError(f'{path}: source coverage bytes differ from sidecar')
            for name in ('levels.csv', 'levels_selection.csv'):
                if metadata['inputs']['sha256'][f'typing/{name}'] != sha(root / name):
                    raise ValueError(f'{path}: source coverage calibration does not match')
            frames.append(CS.indexed(path))
            origins.append({'parquet': CS.file_record(path), 'sidecar': CS.file_record(sidecar),
                            'source_inputs': metadata['inputs']})
        a, b = frames
        CS.equal(a[[AC.INTERIOR_COUNT]], b[[AC.INTERIOR_COUNT]], f'{tag}: source interior raster')
        columns = [f'{prefix}{marker}' for marker in CS.MARKERS for prefix in ('cov__', 'covint__')]
        mixed = b.copy()
        mixed.loc[:, columns] = a.reindex(b.index)[columns]
        AC.write_coverage(run, tag, mixed.reset_index())
        derivations[tag] = {'endocrine': origins[0], 'nonendocrine': origins[1],
                            'selected_endocrine_columns': columns,
                            'output': CS.file_record(run / 'typing/coverage' / f'{tag}.parquet')}
    write_json(run / 'typing/coverage/component_derivation.json', {
        'schema': 'uniform.component_coverage_derivation/1', 'pixels_recomputed': False,
        'method': 'Exact UUID-aligned column selection from hash-verified historical pixel coverage',
        'all_other_columns_from_nonendocrine_source': True, 'images': derivations})
    CC.pin(run, other / 'context', images=tags,
           expected_sha256={tag: sha(other / 'context' / f'{tag}.parquet') for tag in tags},
           source_provenance={'scientific_baseline': source['scientific_baseline'],
                              'selection_reason': args.selection_reason,
                              'package_manifest_sha256': materialized['package_manifest_sha256']})
    CS.pin(run, endocrine, initial, structures, tags, args.selection_reason)
    write_json(run / 'component_staging.json', {
        'schema': 'uniform.component_staging/1', 'run': str(run),
        'classification_control_only': args.classification_control,
        'source_package_manifest_sha256': materialized['package_manifest_sha256'],
        'structure_source': str(structures), 'source_generation_sha256': sha(generation),
        'release_ready': False, 'scientific_accuracy_established': False})
    print(f'Prepared fresh diagnostic component run {run}; identities have not run', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('capture')
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--baseline-description', required=True)
    p.add_argument('--hash-workers', type=int, default=4)
    p.set_defaults(action=capture)
    p = commands.add_parser('materialize')
    p.add_argument('--package', type=Path, required=True)
    p.add_argument('--shared-run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--hash-workers', type=int, default=4)
    p.set_defaults(action=materialize)
    p = commands.add_parser('prepare')
    for name in ('run', 'endocrine-source', 'nonendocrine-source', 'structure-source'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--selection-reason', required=True)
    p.add_argument('--classification-control', action='store_true',
                   help='Typing/audit only: preserve historical initial geometry record without claiming a fresh SG validation')
    p.set_defaults(action=prepare)
    args = parser.parse_args()
    args.action(args)


if __name__ == '__main__':
    main()
