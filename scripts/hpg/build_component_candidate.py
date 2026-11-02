#!/usr/bin/env python3
"""Run an explicit offline component review build inside an allocated HPG job.

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
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ('source_package', 'endocrine_source', 'structure_source', 'donor_source', 'donor_metadata')
OUTPUTS = ('initial_source_run', 'control_run', 'candidate_run', 'assembly_root', 'report_root')


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('CODEX_SITE_SEAL'):
        parser.error('Run through the sealed HPG site wrapper in an allocated SLURM job')
    source_bytes = args.config.read_bytes()
    config = json.loads(source_bytes)
    required = set(INPUTS + OUTPUTS) | {'expected_commit', 'selection_reason'}
    if set(config) != required or not isinstance(config['selection_reason'], str) or not config['selection_reason'].strip():
        parser.error('Config must explicitly name exactly the documented input/output paths, commit and selection reason')
    paths = {name: Path(config[name]).resolve(strict=name in INPUTS) for name in INPUTS + OUTPUTS}
    for name in OUTPUTS:
        path = paths[name]
        if path.exists():
            parser.error(f'{name} already exists; preserve it and select fresh destinations')
        for other_name, other in paths.items():
            if other_name == name:
                continue
            if path == other or path.is_relative_to(other) or other.is_relative_to(path):
                parser.error(f'{name} must be separate from every other input/output: {other_name}')
    def git(*argv):
        return subprocess.check_output(['git', '-C', str(ROOT), *argv], text=True).strip()
    if git('rev-parse', 'HEAD') != config['expected_commit'] or git('status', '--porcelain', '--untracked-files=all'):
        parser.error('Build requires the exact clean, committed source named by the config')
    report = paths['report_root']
    report.mkdir(parents=True)
    (report / 'config.json').write_bytes(source_bytes)
    journal = {'schema': 'uniform.component_build/1', 'status': 'running', 'started_utc': now(),
               'git_commit': config['expected_commit'], 'config_sha256': hashlib.sha256(source_bytes).hexdigest(),
               'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': os.uname().nodename,
               'site_seal': os.environ['CODEX_SITE_SEAL'],
               'site_seal_sha256': os.environ.get('CODEX_SITE_SEAL_SHA256'),
               'release_ready': False, 'scientific_accuracy_established': False, 'steps': []}
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
        step('assemble_app', [py, ROOT / 'scripts/uniform/assemble_uniform.py', '--run', paths['candidate_run'],
             '--donor-source', paths['donor_source'], '--donor-metadata', paths['donor_metadata'],
             '--output', paths['assembly_root'], '--report-root', report / 'assembly',
             '--workers', '6', '--include-collagen-iv'])
        step('verify_exported_states', [py, ROOT / 'scripts/tests/check_uniform_cell_states.py',
             '--run', paths['candidate_run'], '--assembled-root', paths['assembly_root'],
             '--out', report / 'cell_state_exports'])
        if git('rev-parse', 'HEAD') != config['expected_commit'] or git('status', '--porcelain', '--untracked-files=all'):
            raise RuntimeError('Source changed during the build')
        journal.update(status='diagnostic_app_assembled', finished_utc=now())
        save()
        print(f'Offline review build complete: {paths["assembly_root"]}; release acceptance remains blocked.', flush=True)
        return 0
    except Exception as exc:
        journal.update(status='failed', finished_utc=now(), error=f'{type(exc).__name__}: {exc}')
        save()
        raise


if __name__ == '__main__':
    sys.exit(main())
