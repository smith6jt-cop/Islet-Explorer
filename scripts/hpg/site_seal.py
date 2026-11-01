#!/usr/bin/env python3
"""Workstation capture and HiPerGator verify-and-seal for the HPG processing copy (owner, 2026-10-31).

capture (workstation): for every path in a list of absolute paths (regular files and symlinks) record the identity the
typing guards compare (inode, device, size, mtime_ns, ctime_ns), the sha256 taken between two equal stats, and each
symlink's target; with --git-root also git HEAD, refusing a dirty tracked tree. Writes workstation_identities.json,
SHA256SUMS and README.md into a new folder.

verify-and-seal (HPG, inside the site container, so the paths are the workstation's): every file must be a regular
file with the same size and sha256 and the same mtime to the second (Lustre keeps whole seconds); every symlink must
point to the same target; one HPG inode may stand for one workstation identity only. Only then is the seal written,
read-only, for site_identity.py. Any mismatch writes <seal>.failures.json and no seal. Run it without the hook.

    python scripts/hpg/site_seal.py capture --list LIST --out output/hpg_transfer_20261031/capture --git-root .
    python scripts/hpg/site_seal.py verify-and-seal --identities output/hpg_transfer_20261031/capture/workstation_identities.json \
        --identities-sha256 SHA --out output/hpg_transfer_20261031/hpg_site_seal.json
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import socket
import stat as S
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import site_identity as SI  # noqa: E402

IDENTITIES_SCHEMA = 'codex.hpg_workstation_identities/1'


def real_stat(path: str) -> os.stat_result:
    if SI.active():
        raise RuntimeError('site_seal.py must run without the site identity hook (unset CODEX_SITE_SEAL)')
    return os.stat(path)


def identity(st: os.stat_result) -> dict:
    return {'inode': st.st_ino, 'device': st.st_dev, 'mtime_ns': st.st_mtime_ns, 'ctime_ns': st.st_ctime_ns}


def hashed(path: str) -> dict:
    """size, sha256 and identity, refusing a file that changed while it was read."""
    before = real_stat(path)
    if not S.S_ISREG(before.st_mode):
        raise ValueError(f'Not a regular file: {path}')
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 23), b''):
            h.update(block)
    after = real_stat(path)
    if (before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino) != \
            (after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino):
        raise ValueError(f'File changed while hashing: {path}')
    return {'path': path, 'size': after.st_size, 'sha256': h.hexdigest(), 'identity': identity(after)}


def read_list(path: Path) -> list[str]:
    items = [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith('#')]
    bad = [p for p in items if not p.startswith('/') or '\\' in p]
    if bad:
        raise ValueError(f'Paths must be absolute and free of backslashes: {bad[:3]}')
    if len(set(items)) != len(items):
        raise ValueError('Duplicate paths in the list')
    return items


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def capture(a) -> int:
    paths = read_list(Path(a.list))
    files, links = [], []
    for p in paths:
        st = os.lstat(p)
        if S.S_ISLNK(st.st_mode):
            links.append({'path': p, 'target': os.readlink(p)})
        elif S.S_ISREG(st.st_mode):
            files.append(p)
        else:
            raise ValueError(f'Only regular files and symlinks may be listed: {p}')
    git = None
    if a.git_root:
        root = str(Path(a.git_root).resolve())
        dirty = subprocess.run(['git', '-C', root, 'status', '--porcelain', '--untracked-files=no'],
                               capture_output=True, text=True, check=True).stdout
        if dirty.strip():
            raise ValueError('Tracked files are modified; commit before capture:\n' + dirty)
        git = {'root': root, 'head': subprocess.run(['git', '-C', root, 'rev-parse', 'HEAD'], capture_output=True,
                                                    text=True, check=True).stdout.strip()}
    out = Path(a.out)
    out.mkdir(parents=False, exist_ok=False)
    with ProcessPoolExecutor(a.workers) as ex:
        records = sorted(ex.map(hashed, files, chunksize=4), key=lambda r: r['path'])
    doc = {'schema': IDENTITIES_SCHEMA, 'created_utc': datetime.now(timezone.utc).isoformat(),
           'host': socket.gethostname(), 'git': git,
           'list': {'path': str(Path(a.list).resolve()), 'sha256': sha256_file(Path(a.list))},
           'files': records, 'symlinks': sorted(links, key=lambda r: r['path'])}
    target = out / 'workstation_identities.json'
    target.write_text(json.dumps(doc, indent=1) + '\n')
    (out / 'SHA256SUMS').write_text(''.join(f"{r['sha256']}  {r['path']}\n" for r in records))
    total = sum(r['size'] for r in records)
    (out / 'README.md').write_text(
        f'# HPG transfer capture\n\nWritten by `scripts/hpg/site_seal.py capture` on {doc["host"]} at '
        f'{doc["created_utc"]}: {len(records)} files ({total / 1e9:.1f} GB) and {len(links)} symlinks from '
        f'`{doc["list"]["path"]}`, git HEAD {git and git["head"]}. `SHA256SUMS` lists every file; '
        '`workstation_identities.json` also holds the identities the typing guards compare. On HPG, '
        '`site_seal.py verify-and-seal` checks the copy against it and writes the seal used by '
        '`scripts/hpg/site_identity.py`. See docs/hpg_processing_20261031.md.\n')
    print(f'{len(records)} files, {total / 1e9:.1f} GB, {len(links)} symlinks')
    print(f'workstation_identities.json sha256 {sha256_file(target)}')
    return 0


def verify_and_seal(a) -> int:
    source = Path(a.identities)
    actual = sha256_file(source)
    if actual != a.identities_sha256:
        raise ValueError(f'{source} has sha256 {actual}, expected {a.identities_sha256}')
    doc = json.loads(source.read_text())
    if doc.get('schema') != IDENTITIES_SCHEMA:
        raise ValueError(f'Unsupported identities schema: {doc.get("schema")}')
    out = Path(a.out)
    if out.exists():
        raise FileExistsError(out)
    failures = []
    for link in doc['symlinks']:
        p = link['path']
        if not os.path.islink(p):
            failures.append({'path': p, 'problem': 'symlink missing'})
        elif os.readlink(p) != link['target']:
            failures.append({'path': p, 'problem': 'symlink target differs', 'target': os.readlink(p)})
    expected = {r['path']: r for r in doc['files']}
    present = [p for p in expected if os.path.isfile(p) and not os.path.islink(p)]
    absent = set(expected) - set(present)
    failures += [{'path': p, 'problem': 'file missing or not a regular file'} for p in sorted(absent)]
    with ProcessPoolExecutor(a.workers) as ex:
        found = list(ex.map(hashed, present, chunksize=4))
    files, by_inode = [], {}
    for r in found:
        w = expected[r['path']]
        problems = [k for k in ('size', 'sha256') if r[k] != w[k]]
        if r['identity']['mtime_ns'] // 10 ** 9 != w['identity']['mtime_ns'] // 10 ** 9:
            problems.append('mtime seconds')
        if problems:
            failures.append({'path': r['path'], 'problem': ', '.join(problems) + ' differ'})
            continue
        h, ws = r['identity'], w['identity']
        if by_inode.setdefault(h['inode'], ws) != ws:
            failures.append({'path': r['path'], 'problem': 'HPG inode already stands for another workstation identity'})
            continue
        files.append({'path': r['path'], 'size': r['size'], 'sha256': r['sha256'], 'hpg': h, 'workstation': ws})
    if failures:
        report = out.with_name(out.name + '.failures.json')
        report.write_text(json.dumps({'identities': str(source), 'failures': failures}, indent=1) + '\n')
        print(f'{len(failures)} failures; no seal written; see {report}', file=sys.stderr)
        return 2
    seal = {'schema': SI.SCHEMA, 'created_utc': datetime.now(timezone.utc).isoformat(), 'host': socket.gethostname(),
            'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'container': os.environ.get('APPTAINER_CONTAINER') or os.environ.get('SINGULARITY_CONTAINER'),
            'workstation_identities': {'path': str(source.resolve()), 'sha256': actual, 'git': doc['git']},
            'symlinks_checked': len(doc['symlinks']), 'files': sorted(files, key=lambda f: f['path'])}
    with out.open('x') as stream:
        stream.write(json.dumps(seal, indent=1) + '\n')
    out.chmod(0o444)
    print(f'{len(files)} files and {len(doc["symlinks"])} symlinks verified')
    print(f'seal {out} sha256 {sha256_file(out)}')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='command', required=True)
    c = sub.add_parser('capture')
    c.add_argument('--list', required=True)
    c.add_argument('--out', required=True)
    c.add_argument('--git-root')
    c.add_argument('--workers', type=int, default=4)
    v = sub.add_parser('verify-and-seal')
    v.add_argument('--identities', required=True)
    v.add_argument('--identities-sha256', required=True)
    v.add_argument('--out', required=True)
    v.add_argument('--workers', type=int, default=16)
    a = ap.parse_args()
    return capture(a) if a.command == 'capture' else verify_and_seal(a)


if __name__ == '__main__':
    raise SystemExit(main())
