#!/usr/bin/env python3
"""Create a fresh blinded native QuPath sample with a separate prediction key.

The final identity table is the sampling population (including any final Excluded stratum,
but not QC-ineligible cells absent from identity). Previously reviewed cells are excluded.
No accuracy estimate is computed: reviewer labels must be supplied independently first.

    python scripts/uniform/blind_review_project.py --run <candidate-run> \\
      --project <new-blind-project> --key-dir <separate-private-key-directory>
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE / 'typing'))
import codex_typing as CT

QUPATH = '/home/smith6jt/QuPath/bin/QuPath'
SCRIPT = HERE / 'qupath' / 'import_blind_review.groovy'
DEFAULT_EXCLUSIONS = REPO / 'output' / 'typing_accuracy_20261029' / 'manifest.csv'
PUBLIC_COLUMNS = ['sample_id', 'image', 'x_um', 'y_um', 'reviewer_label', 'reviewer_confidence', 'notes']


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def exclusions(path: Path) -> set[tuple[str, str]]:
    table = pd.read_parquet(path) if path.suffix == '.parquet' else pd.read_csv(path, dtype=str)
    if not {'image', 'object_id'} <= set(table):
        raise ValueError('Exclusion manifest must contain image and object_id')
    if table[['image', 'object_id']].isna().any().any() or table[['image', 'object_id']].eq('').any().any():
        raise ValueError('Exclusion manifest contains missing image/object IDs')
    return set(table[['image', 'object_id']].itertuples(index=False, name=None))


def sample_by_label(population: pd.DataFrame, excluded: set[tuple[str, str]],
                    per_label: int = 20, seed: int = 1) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stratified without replacement; randomized presentation order never encodes the label."""
    if per_label < 1:
        raise ValueError('per_label must be positive')
    columns = ['image', 'object_id', 'label']
    frame = population[columns].copy()
    if frame.empty:
        raise ValueError('Sampling population is empty')
    if frame.isna().any().any() or frame.eq('').any().any():
        raise ValueError('Sampling population has missing IDs/labels')
    if frame[['image', 'object_id']].duplicated().any():
        raise ValueError('Sampling population has duplicate image/object IDs')
    frame = frame.sort_values(['image', 'object_id']).reset_index(drop=True)
    blocked = np.array([(image, oid) in excluded for image, oid in frame[['image', 'object_id']].itertuples(index=False)])
    eligible = frame.loc[~blocked]
    if eligible.empty:
        raise ValueError('No previously unreviewed cells remain')
    generator = np.random.default_rng(seed)
    parts, strata = [], []
    for label, group in frame.groupby('label', sort=True):
        available = eligible[eligible['label'] == label]
        count = min(per_label, len(available))
        info = {'label': label, 'full_population_n': len(group), 'excluded_previous_n': len(group) - len(available),
                'eligible_population_n': len(available), 'sample_n': count,
                'full_prevalence': len(group) / len(frame), 'eligible_prevalence': len(available) / len(eligible),
                'inclusion_probability': count / len(available) if len(available) else 0.0,
                'inverse_sampling_weight': len(available) / count if count else np.nan}
        strata.append(info)
        if count:
            sampled = available.iloc[generator.choice(len(available), count, replace=False)].copy()
            for key, value in info.items():
                if key != 'label':
                    sampled[key] = value
            parts.append(sampled)
    sampled = pd.concat(parts, ignore_index=True)
    sampled = sampled.iloc[generator.permutation(len(sampled))].reset_index(drop=True)
    sampled.insert(0, 'sample_id', [f'Review-{i + 1:06d}' for i in range(len(sampled))])
    sampled['blind_uuid'] = [str(uuid.uuid5(uuid.NAMESPACE_URL, f'blind-review/{seed}/{row.image}/{row.object_id}'))
                            for row in sampled.itertuples()]
    return sampled, pd.DataFrame(strata)


def public_table(key: pd.DataFrame) -> pd.DataFrame:
    table = key[['sample_id', 'review_image', 'x_um', 'y_um']].rename(columns={'review_image': 'image'}).copy()
    table['reviewer_label'] = ''
    table['reviewer_confidence'] = ''
    table['notes'] = ''
    return table[PUBLIC_COLUMNS]


def minimal_project(source: dict, destination: Path, image_names: dict[str, str]) -> dict:
    """Keep only raw server definitions and selected entries; do not copy project metadata/classifiers/scripts."""
    images = []
    for entry in source['images']:
        if entry['imageName'] in image_names:
            images.append({'entryID': entry['entryID'], 'imageName': image_names[entry['imageName']],
                           'serverBuilder': copy.deepcopy(entry['serverBuilder'])})
    if len(images) != len(image_names):
        raise ValueError('Some sampled images are missing from the source project')
    # QuPath 0.7 requires both scalar timestamps while loading the project. They
    # carry no prediction metadata; omitting them prevents the sanitizer running.
    return {'version': source['version'], 'createTimestamp': source['createTimestamp'],
            'modifyTimestamp': source['modifyTimestamp'], 'uri': destination.as_uri(), 'images': images,
            'lastID': max(int(e['entryID']) for e in images)}


def ensure_separate(project: Path, key_dir: Path, run: Path) -> None:
    if project == key_dir or project in key_dir.parents or key_dir in project.parents:
        raise ValueError('Blind project and private key directory must be separate sibling locations')
    if any(output == run or output in run.parents or run in output.parents for output in (project, key_dir)):
        raise ValueError('Output paths cannot contain, replace, or be inside the source run')
    if project.exists() or key_dir.exists():
        raise FileExistsError('Use fresh project and key directories; existing review data is never overwritten')


def build(run: Path, project: Path, key_dir: Path, exclude_path: Path,
          per_label: int = 20, seed: int = 1, workers: int = 2) -> dict:
    ensure_separate(project, key_dir, run)
    if workers < 1:
        raise ValueError('workers must be positive')
    rows = CT.cohort_rows()
    tables, inputs = [], {str(exclude_path): sha256(exclude_path)}
    for row in rows:
        tag = row['image_name']
        path = run / 'typing' / 'identity' / f'{tag}.parquet'
        table = pd.read_parquet(path, columns=['object_id', 'label'])
        table.insert(0, 'image', tag)
        tables.append(table)
        inputs[str(path)] = sha256(path)
    excluded = exclusions(exclude_path)
    key, strata = sample_by_label(pd.concat(tables, ignore_index=True), excluded, per_label, seed)
    image_names = {tag: f'Review image {i + 1:03d}' for i, tag in enumerate(sorted(key['image'].unique()))}
    key['review_image'] = key['image'].map(image_names)
    locations = []
    for tag in image_names:
        path = run / 'cells' / f'{tag}.parquet'
        table = pd.read_parquet(path, columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
        table = table.rename(columns={'Object ID': 'object_id', 'Centroid X µm': 'x_um', 'Centroid Y µm': 'y_um'})
        table['image'] = tag
        locations.append(table)
        inputs[str(path)] = sha256(path)
    key = key.merge(pd.concat(locations, ignore_index=True), on=['image', 'object_id'], how='left', validate='one_to_one')
    if key[['x_um', 'y_um']].isna().any().any():
        raise ValueError('Sampled cells lack native image coordinates')
    source_path = run / 'qupath' / 'project.qpproj'
    source = json.loads(source_path.read_text())
    inputs[str(source_path)] = sha256(source_path)
    key_dir.mkdir(parents=True, mode=0o700)
    staging = key_dir / '_project_in_progress'
    staging.mkdir()
    entries = minimal_project(source, staging / 'project.qpproj', image_names)
    (staging / 'project.qpproj').write_text(json.dumps(entries, indent=2) + '\n')
    private_assignments = key_dir / 'assignments'
    private_assignments.mkdir()
    logs = key_dir / 'logs'; logs.mkdir()
    source_entries = {e['imageName']: e for e in source['images']}
    for tag, neutral in image_names.items():
        entry_id = str(source_entries[tag]['entryID'])
        source_data = run / 'qupath' / 'data' / entry_id / 'data.qpdata'
        inputs[str(source_data)] = sha256(source_data)
        destination = staging / 'data' / entry_id
        destination.mkdir(parents=True)
        shutil.copy2(source_data, destination / 'data.qpdata')
        subset = key[key['image'] == tag]
        assignments = subset[['object_id', 'sample_id', 'blind_uuid']].to_dict('records')
        (private_assignments / f'{neutral}.json').write_text(json.dumps(assignments, indent=2) + '\n')
    key.to_csv(key_dir / 'prediction_key.csv', index=False)
    strata.to_csv(key_dir / 'sampling_strata.csv', index=False)
    instructions = ('Review the named cells against their raw image channels. Enter an independent label for each sample_id '
                    'in reviewer_labels.csv; leave an uncertain cell unresolved and explain in notes. All sample cells have '
                    'the same neutral class. Original cell/nucleus outlines and raw images are preserved. Cell measurements, '
                    'prediction metadata and endocrine/structure annotations are hidden.\n'
                    'No accuracy estimate exists until independent reviewer labels are supplied.\n')
    public_table(key).to_csv(staging / 'reviewer_labels.csv', index=False)
    (staging / 'REVIEW.txt').write_text(instructions)
    def import_one(neutral):
        log = logs / f'{neutral}.log'
        assignment = private_assignments / f'{neutral}.json'
        report = key_dir / f'{neutral}.sanitization.json'
        with log.open('w') as stream:
            result = subprocess.run([QUPATH, 'script', '-p', str(staging / 'project.qpproj'), '-i', neutral,
                                     '-a', str(assignment), '-a', str(report), str(SCRIPT)],
                                    stdout=stream, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                    env={**os.environ, 'JAVA_TOOL_OPTIONS': '-Xmx64g -XX:ActiveProcessorCount=12'})
        if result.returncode or not report.is_file():
            raise RuntimeError(f'Blind import failed for {neutral}: {log}')
        data = json.loads(report.read_text())
        if not data.get('passed'):
            raise ValueError(f'Blind sanitization failed for {neutral}: {report}')
        return data
    with ThreadPoolExecutor(workers) as pool:
        checks = list(pool.map(import_one, image_names.values()))
    # A QuPath save creates a backup of the original, unblinded data. Never deliver those backups.
    for entry in entries['images']:
        folder = staging / 'data' / str(entry['entryID'])
        for path in folder.iterdir():
            if path.name != 'data.qpdata':
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
    (staging / 'classifiers').mkdir(exist_ok=True)
    (staging / 'classifiers' / 'classes.json').write_text(json.dumps({'pathClasses': [
        {'name': 'Tissue', 'color': -10053223}, {'name': 'Review cell', 'color': -8355712}]}, indent=2) + '\n')
    for path, digest in inputs.items():
        if sha256(Path(path)) != digest:
            raise ValueError(f'Source changed during blind project creation: {path}')
    entries['uri'] = (project / 'project.qpproj').as_uri()
    (staging / 'project.qpproj').write_text(json.dumps(entries, indent=2) + '\n')
    public_manifest = {'status': 'awaiting_independent_manual_labels', 'samples': len(key), 'images': len(image_names),
                       'accuracy_computed': False, 'seed': seed, 'per_label_target': per_label,
                       'blinding': 'neutral cell classes, randomized IDs/order; no predictions or private key in project'}
    (staging / 'blind_review_manifest.json').write_text(json.dumps(public_manifest, indent=2) + '\n')
    private_manifest = {**public_manifest, 'source_run': str(run), 'project': str(project),
                        'sampling_population': 'final identity rows excluding prior review image/object IDs',
                        'weight_estimand': 'eligible previously unreviewed population; inverse inclusion probability N_h/n_h',
                        'excluded_manifest': str(exclude_path), 'excluded_manifest_sha256': sha256(exclude_path),
                        'source_inputs': inputs, 'sanitization': checks, 'script_sha256': sha256(SCRIPT),
                        'driver_sha256': sha256(Path(__file__))}
    (key_dir / 'manifest.json').write_text(json.dumps(private_manifest, indent=2) + '\n')
    project.parent.mkdir(parents=True, exist_ok=True)
    if project.exists():
        raise FileExistsError(f'Blind project appeared during creation; refusing to overwrite: {project}')
    staging.rename(project)
    return private_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--key-dir', type=Path, required=True)
    parser.add_argument('--exclude-manifest', type=Path, default=DEFAULT_EXCLUSIONS)
    parser.add_argument('--per-label', type=int, default=20)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    report = build(args.run.resolve(), args.project.resolve(), args.key_dir.resolve(), args.exclude_manifest.resolve(),
                   args.per_label, args.seed, args.workers)
    print(json.dumps({k: report[k] for k in ('status', 'samples', 'images', 'project', 'accuracy_computed')}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
