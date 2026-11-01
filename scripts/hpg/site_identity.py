"""HiPerGator site identity: sealed copies answer os.stat with their workstation file identity (owner, 2026-10-31).

The typing guards record a file's {path, size, mtime_ns, ctime_ns, inode} and compare it exactly
(structure_generation._stat, new_detection_filter.file_stat, refresh_ready_cells.file_stat,
run_remaining_qc_pipeline.stat_identity, final_reconciliation, collagen_structures' raw image). No copy keeps ctime or
inode, and HPG's Lustre /blue keeps mtime to the second. Four of those functions sit in files whose own sha256 is
pinned in records, so none of them is edited. Instead, on HPG only, sitecustomize.py wraps os.stat (pathlib's
Path.stat calls it): a regular file whose current (inode, size, mtime_ns, ctime_ns) equals a sealed tuple is reported
with its workstation inode, device and times. Everything else is reported as is: a changed, replaced or new file, a
directory, lstat, fstat, DirEntry.stat and dir_fd calls. A change after sealing therefore still fails the guards, and
a process started without the hook (python -S, -I or -E) fails them too. The seal is written by
site_seal.py verify-and-seal after the copy's sha256, size and mtime (to the second) match the workstation's, and is
pinned by CODEX_SITE_SEAL_SHA256. See docs/hpg_processing_20261031.md.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import socket
import stat as S
import sys

SCHEMA = 'codex.hpg_site_seal/1'
EXIT_CODE = 78
HERE = Path(__file__).resolve().parent


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 23), b''):
            h.update(block)
    return h.hexdigest()


def load_table(seal_path, expected_sha256: str) -> tuple[dict, dict]:
    """{HPG inode: (size, HPG mtime_ns, HPG ctime_ns, workstation inode, device, mtime_ns, ctime_ns)} and the seal's
    summary. Refuses a seal whose bytes do not hash to expected_sha256."""
    raw = Path(seal_path).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected_sha256:
        raise ValueError(f'Site seal {seal_path} has sha256 {actual}, expected {expected_sha256}')
    seal = json.loads(raw)
    if seal.get('schema') != SCHEMA:
        raise ValueError(f'Unsupported site seal schema: {seal.get("schema")}')
    table = {}
    for f in seal['files']:
        h, w = f['hpg'], f['workstation']
        value = (f['size'], h['mtime_ns'], h['ctime_ns'], w['inode'], w['device'], w['mtime_ns'], w['ctime_ns'])
        if table.setdefault(h['inode'], value) != value:
            raise ValueError(f'Site seal maps HPG inode {h["inode"]} to two workstation identities')
    info = {'seal': str(Path(seal_path).resolve()), 'seal_sha256': actual, 'files': len(seal['files']),
            'workstation_identities_sha256': seal['workstation_identities']['sha256'],
            'hook_sha256': {p.name: sha256(p) for p in (HERE / 'sitecustomize.py', HERE / 'site_identity.py')}}
    return table, info


def translated(st: os.stat_result, entry: tuple) -> os.stat_result:
    """st with the workstation inode, device, mtime and ctime of entry; every other field unchanged."""
    cls, (seq, extra) = st.__reduce__()
    seq, extra = list(seq), dict(extra)
    _, _, _, ino, dev, mtime_ns, ctime_ns = entry
    seq[1], seq[2] = ino, dev
    for i, name, ns in ((8, 'st_mtime', mtime_ns), (9, 'st_ctime', ctime_ns)):
        sec, nsec = divmod(ns, 10 ** 9)
        seq[i] = sec
        extra[name] = sec + nsec * 1e-9  # as CPython's fill_time
        extra[name + '_ns'] = ns
    return cls(seq, extra)


def wrap(real_stat, table: dict, info: dict):
    def stat(path, *, dir_fd=None, follow_symlinks=True):
        st = real_stat(path, dir_fd=dir_fd, follow_symlinks=follow_symlinks)
        if dir_fd is None and follow_symlinks and not isinstance(path, int) and S.S_ISREG(st.st_mode):
            entry = table.get(st.st_ino)
            if entry is not None and entry[:3] == (st.st_size, st.st_mtime_ns, st.st_ctime_ns):
                return translated(st, entry)
        return st
    stat.codex_site = info
    stat.__wrapped__ = real_stat
    stat.__doc__ = real_stat.__doc__
    return stat


def install_from_environment() -> bool:
    """Called by sitecustomize.py; any failure ends the process with EXIT_CODE (site.py would swallow exceptions)."""
    path = os.environ.get('CODEX_SITE_SEAL')
    if not path:
        return False
    try:
        if hasattr(os.stat, 'codex_site'):
            return True
        table, info = load_table(path, os.environ.get('CODEX_SITE_SEAL_SHA256', ''))
        os.stat = wrap(os.stat, table, info)
        return True
    except BaseException as e:
        print(f'CODEX site identity hook refused: {type(e).__name__}: {e}', file=sys.stderr, flush=True)
        os._exit(EXIT_CODE)


def active() -> dict | None:
    return getattr(os.stat, 'codex_site', None)


def journal_entry() -> dict | None:
    """The site record for a driver's journal: None on the workstation; refuses under SLURM without the hook."""
    info = active()
    if info is None:
        if os.environ.get('SLURM_JOB_ID'):
            raise RuntimeError('Under SLURM the HPG site identity hook must be active (CODEX_SITE_SEAL, PYTHONPATH)')
        return None
    return {**info, 'host': socket.gethostname(), 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'container': os.environ.get('APPTAINER_CONTAINER') or os.environ.get('SINGULARITY_CONTAINER')}
