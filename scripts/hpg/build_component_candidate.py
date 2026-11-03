#!/usr/bin/env python3
"""Build and check an explicit component candidate inside an allocated HPG job.

Every scientific destination must be fresh. Failed output is retained, and no
release attestation, deployment, service change or automatic retry is performed.
The immutable config names every input and output; there are no version defaults.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ('source_package', 'endocrine_source', 'structure_source', 'donor_source', 'donor_metadata',
          'classifier_run', 'classifier_replay', 'geometry_reference', 'geometry_provenance')
OUTPUTS = ('initial_source_run', 'control_run', 'candidate_run', 'assembly_root', 'report_root')
PINNED_REPORTS = ('classifier_replay', 'geometry_provenance')


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def validate_config(config):
    """Check all destinations and explicit receipt pins before creating output."""
    required = set(INPUTS + OUTPUTS) | {'expected_commit', 'selection_reason'} | {
        name + '_sha256' for name in PINNED_REPORTS}
    if set(config) != required:
        raise ValueError('Config must explicitly name exactly the documented paths, commit, selection reason and receipt hashes')
    if not isinstance(config['selection_reason'], str) or not config['selection_reason'].strip():
        raise ValueError('A nonempty selection reason is required')
    if not isinstance(config['expected_commit'], str) or not re.fullmatch(r'[0-9a-f]{40}', config['expected_commit']):
        raise ValueError('expected_commit must be a full Git commit SHA')
    if any(not isinstance(config[name], str) or not Path(config[name]).is_absolute() for name in INPUTS + OUTPUTS):
        raise ValueError('Every input and output must be an explicit absolute path')
    # Check the unresolved pathname too: resolve() would hide a dangling link.
    for name in OUTPUTS:
        if Path(config[name]).exists() or Path(config[name]).is_symlink():
            raise ValueError(f'{name} already exists; preserve it and select fresh destinations')
    paths = {name: Path(config[name]).resolve(strict=name in INPUTS) for name in INPUTS + OUTPUTS}
    for name in OUTPUTS:
        for other_name, other in paths.items():
            if other_name == name:
                continue
            path = paths[name]
            if path == other or path.is_relative_to(other) or other.is_relative_to(path):
                raise ValueError(f'{name} must be separate from every other input/output: {other_name}')
    for name in PINNED_REPORTS:
        digest = config[name + '_sha256']
        if not isinstance(digest, str) or not re.fullmatch(r'[0-9a-f]{64}', digest):
            raise ValueError(f'{name}_sha256 must explicitly pin the receipt bytes')
        if sha256(paths[name]) != digest:
            raise ValueError(f'{name}: receipt SHA256 differs from the explicit config pin')
    return paths


def checked_gate(path, *, schema=None, pass_key, pass_value=True, expected=None):
    """A zero process exit is insufficient: bind the actual passing report."""
    data = json.loads(path.read_bytes())
    if schema is not None and data.get('schema') != schema:
        raise ValueError(f'{path}: unexpected gate schema')
    if type(data.get(pass_key)) is not type(pass_value) or data.get(pass_key) != pass_value:
        raise ValueError(f'{path}: gate did not pass')
    for key, value in (expected or {}).items():
        if data.get(key) != str(value):
            raise ValueError(f'{path}: gate {key} differs from this build')
    images = data.get('images')
    if not isinstance(images, (list, dict)) or len(images) != 20:
        raise ValueError(f'{path}: gate must cover all 20 cohort images')
    if any(not isinstance(image, str) or not image for image in images) or len(set(images)) != 20:
        raise ValueError(f'{path}: gate has invalid or duplicate image IDs')
    return {'path': str(path), 'sha256': sha256(path), 'passed': True, 'images': sorted(images)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('CODEX_SITE_SEAL'):
        parser.error('Run through the sealed HPG site wrapper in an allocated SLURM job')
    source_bytes = args.config.read_bytes()
    config = json.loads(source_bytes)
    try:
        paths = validate_config(config)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    def git(*argv):
        return subprocess.check_output(['git', '-C', str(ROOT), *argv], text=True).strip()
    if git('rev-parse', 'HEAD') != config['expected_commit'] or git('status', '--porcelain', '--untracked-files=all'):
        parser.error('Build requires the exact clean, committed source named by the config')
    report = paths['report_root']
    report.mkdir(parents=True)
    (report / 'config.json').write_bytes(source_bytes)
    journal = {'schema': 'uniform.component_build/2', 'status': 'running', 'started_utc': now(),
               'git_commit': config['expected_commit'], 'config_sha256': hashlib.sha256(source_bytes).hexdigest(),
               'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': os.uname().nodename,
               'site_seal': os.environ['CODEX_SITE_SEAL'],
               'site_seal_sha256': os.environ.get('CODEX_SITE_SEAL_SHA256'),
               'release_ready': False, 'scientific_accuracy_established': False, 'gates': {}, 'steps': []}
    def save():
        target = report / 'build_journal.json'
        temporary = target.with_suffix('.tmp')
        temporary.write_text(json.dumps(journal, indent=2) + '\n')
        temporary.replace(target)
    def step(name, argv):
        entry = {'name': name, 'command': [str(x) for x in argv], 'started_utc': now()}
        journal['steps'].append(entry)
        save()
        print(f'[{entry["started_utc"]}] {name}', flush=True)
        result = subprocess.run(entry['command'], cwd=ROOT, check=False)
        entry.update(exit_code=result.returncode, finished_utc=now())
        save()
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, entry['command'])
    py = sys.executable
    stage = ROOT / 'scripts/hpg/stage_component_build.py'
    try:
        step('materialize_initial_source', [py, stage, 'materialize', '--package', paths['source_package'],
             '--shared-run', paths['structure_source'], '--output', paths['initial_source_run'], '--hash-workers', '4'])
        for name, structure, extra in (
            ('control_run', paths['initial_source_run'], ['--classification-control']),
            ('candidate_run', paths['structure_source'], []),
        ):
            step('prepare_' + name, [py, stage, 'prepare', '--run', paths[name],
                 '--endocrine-source', paths['endocrine_source'], '--nonendocrine-source', paths['initial_source_run'],
                 '--structure-source', structure, '--selection-reason', config['selection_reason'], *extra])
            step('identity_' + name, [py, ROOT / 'scripts/uniform/typing/identity.py',
                 '--run', paths[name], '--workers', '4'])
        step('lymphatic_invariance', [py, ROOT / 'scripts/uniform/audit_structure_update.py',
             '--baseline', paths['control_run'], '--candidate', paths['candidate_run'],
             '--output', report / 'lymphatic_invariance', '--images', 'all', '--workers', '4'])
        journal['gates']['lymphatic_invariance'] = checked_gate(report / 'lymphatic_invariance/summary.json',
             schema='uniform.structure_update_invariant_audit/1', pass_key='passed',
             expected={'baseline': paths['control_run'], 'candidate': paths['candidate_run']})
        save()
        object_report = report / 'endocrine_object_audit/summary.json'
        step('endocrine_classifier_lineage', [py, ROOT / 'scripts/uniform/typing/endocrine_object_audit.py',
             '--run', paths['candidate_run'], '--out', object_report.parent,
             '--classifier-run', paths['classifier_run'], '--classifier-replay', paths['classifier_replay'],
             '--classifier-replay-sha256', config['classifier_replay_sha256'],
             '--geometry-reference', paths['geometry_reference'], '--geometry-provenance', paths['geometry_provenance'],
             '--geometry-provenance-sha256', config['geometry_provenance_sha256']])
        journal['gates']['endocrine_classifier_lineage'] = checked_gate(object_report,
             schema='uniform.endocrine_object_audit/2', pass_key='passes',
             expected={'audit_mode': 'classifier_lineage', 'run': paths['candidate_run']})
        save()
        step('assemble_app', [py, ROOT / 'scripts/uniform/assemble_uniform.py', '--run', paths['candidate_run'],
             '--donor-source', paths['donor_source'], '--donor-metadata', paths['donor_metadata'],
             '--output', paths['assembly_root'], '--report-root', report / 'assembly',
             '--endocrine-object-audit', object_report, '--workers', '6', '--include-collagen-iv'])
        step('verify_exported_states', [py, ROOT / 'scripts/tests/check_uniform_cell_states.py',
             '--run', paths['candidate_run'], '--assembled-root', paths['assembly_root'],
             '--out', report / 'cell_state_exports'])
        journal['gates']['exported_cell_states'] = checked_gate(report / 'cell_state_exports/summary.json',
             pass_key='result', pass_value='pass',
             expected={'run': paths['candidate_run'], 'assembled_root': paths['assembly_root']})
        if len({tuple(gate['images']) for gate in journal['gates'].values()}) != 1:
            raise ValueError('The build gates cover different image inventories')
        if git('rev-parse', 'HEAD') != config['expected_commit'] or git('status', '--porcelain', '--untracked-files=all'):
            raise RuntimeError('Source changed during the build')
        for name in PINNED_REPORTS:
            if sha256(paths[name]) != config[name + '_sha256']:
                raise RuntimeError(f'{name}: pinned source receipt changed during the build')
        journal.update(status='app_assembled', finished_utc=now())
        save()
        print(f'App assembled and build gates passed: {paths["assembly_root"]}; release and host checks remain separate.', flush=True)
        return 0
    except Exception as exc:
        journal.update(status='failed', finished_utc=now(), error=f'{type(exc).__name__}: {exc}')
        save()
        raise


if __name__ == '__main__':
    sys.exit(main())
