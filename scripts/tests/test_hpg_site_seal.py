"""The HPG site identity hook answers os.stat with the workstation identity only for sealed, unchanged files."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

HPG = Path(__file__).resolve().parents[1] / 'hpg'
UNIFORM = Path(__file__).resolve().parents[1] / 'uniform'
sys.path.insert(0, str(HPG))
sys.path.insert(0, str(UNIFORM))
sys.path.insert(0, str(UNIFORM / 'typing'))
import site_identity as SI  # noqa: E402

WS_DEVICE = 2049
REAL_STAT = os.stat


def ws_identity(st, shift):
    return {'inode': st.st_ino + shift, 'device': WS_DEVICE, 'mtime_ns': st.st_mtime_ns // 10 ** 9 * 10 ** 9 + 840296517,
            'ctime_ns': st.st_ctime_ns - 10 ** 12 + 123}


def make_seal(tmp_path, paths, shift=10 ** 6):
    entries = []
    for p in paths:
        st = REAL_STAT(p)
        entries.append({'path': str(p), 'size': st.st_size, 'sha256': SI.sha256(p),
                        'hpg': {'inode': st.st_ino, 'device': st.st_dev, 'mtime_ns': st.st_mtime_ns,
                                'ctime_ns': st.st_ctime_ns},
                        'workstation': ws_identity(st, shift)})
    seal = tmp_path / 'seal.json'
    seal.write_text(json.dumps({'schema': SI.SCHEMA, 'workstation_identities': {'sha256': '0' * 64}, 'files': entries}))
    return seal, SI.sha256(seal)


@pytest.fixture
def sealed(tmp_path, monkeypatch):
    f = tmp_path / 'data' / 'a.bin'
    f.parent.mkdir()
    f.write_bytes(b'abc' * 1000)
    seal, digest = make_seal(tmp_path, [f])
    table, info = SI.load_table(seal, digest)
    monkeypatch.setattr(os, 'stat', SI.wrap(REAL_STAT, table, info))
    time.sleep(0.02)  # kernel timestamps advance per clock tick; let a later change be visible
    return f, json.loads(seal.read_text())['files'][0]['workstation']


def run_python(code, env_extra, *args, flags=()):
    env = {k: v for k, v in os.environ.items() if not k.startswith('CODEX_SITE')}
    env.update(PYTHONPATH=str(HPG), **env_extra)
    return subprocess.run([sys.executable, *flags, '-c', code, *args], capture_output=True, text=True, env=env)


def test_inert_without_seal_variable():
    r = run_python('import os, posix; print(os.stat is posix.stat)', {})
    assert r.returncode == 0 and r.stdout.strip() == 'True'


def test_wrong_seal_sha_ends_the_process(tmp_path):
    f = tmp_path / 'a'
    f.write_text('x')
    seal, _ = make_seal(tmp_path, [f])
    r = run_python('print("ran")', {'CODEX_SITE_SEAL': str(seal), 'CODEX_SITE_SEAL_SHA256': 'f' * 64})
    assert r.returncode == SI.EXIT_CODE and 'ran' not in r.stdout


def test_sealed_file_reports_workstation_identity(sealed):
    f, ws = sealed
    st = f.stat()
    real = REAL_STAT(f)
    assert (st.st_ino, st.st_dev, st.st_mtime_ns, st.st_ctime_ns) == \
           (ws['inode'], ws['device'], ws['mtime_ns'], ws['ctime_ns'])
    assert (st.st_size, st.st_mode, st.st_nlink, st.st_uid) == (real.st_size, real.st_mode, real.st_nlink, real.st_uid)
    assert st.st_mtime == ws['mtime_ns'] // 10 ** 9 + ws['mtime_ns'] % 10 ** 9 * 1e-9
    assert st[8] == ws['mtime_ns'] // 10 ** 9 and st[9] == ws['ctime_ns'] // 10 ** 9


@pytest.mark.parametrize('change', ['write', 'utime', 'rename', 'hardlink', 'chmod'])
def test_any_change_reports_the_real_identity(sealed, change):
    f, ws = sealed
    if change == 'write':
        f.write_bytes(b'abc' * 1000)
    elif change == 'utime':
        os.utime(f)
    elif change == 'rename':
        f = f.rename(f.with_name('b.bin'))
    elif change == 'hardlink':
        os.link(f, f.with_name('link.bin'))
    else:
        f.chmod(0o600)
    assert os.stat(f).st_ino != ws['inode']
    assert f.stat() == REAL_STAT(f)


def test_directories_lstat_and_new_files_are_untranslated(sealed):
    f, ws = sealed
    new = f.with_name('new.bin')
    new.write_bytes(b'abc' * 1000)
    assert os.stat(new) == REAL_STAT(new)
    assert os.stat(f.parent) == REAL_STAT(f.parent)
    assert os.lstat(f).st_ino != ws['inode']
    assert os.stat(f, follow_symlinks=False).st_ino != ws['inode']


def test_symlink_followed_to_a_sealed_file_is_translated(sealed):
    f, ws = sealed
    link = f.with_name('alias')
    link.symlink_to(f)
    assert link.stat().st_ino == ws['inode']


def test_hard_links_share_one_identity_and_conflicts_are_refused(tmp_path):
    f = tmp_path / 'a'
    f.write_text('x')
    g = tmp_path / 'b'
    os.link(f, g)
    seal, digest = make_seal(tmp_path, [f, g])
    table, _ = SI.load_table(seal, digest)
    assert len(table) == 1
    doc = json.loads(seal.read_text())
    doc['files'][1]['workstation']['inode'] += 1
    seal.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match='two workstation identities'):
        SI.load_table(seal, SI.sha256(seal))


def test_samestat_never_matches_translated_and_real(sealed):
    f, _ = sealed
    other = f.with_name('other.bin')
    other.write_text('y')
    assert os.path.samefile(f, f)
    assert not os.path.samestat(os.stat(f), REAL_STAT(other))
    with open(f, 'rb') as stream:
        assert not os.path.samestat(os.stat(f), os.fstat(stream.fileno()))


def test_guard_identity_functions_see_the_workstation_identity(sealed):
    import new_detection_filter as NF
    import refresh_ready_cells as RF
    import structure_generation as SG
    f, ws = sealed
    for record in (SG._stat(f), NF.file_stat(f), RF.file_stat(f)):
        assert (record['inode'], record['mtime_ns'], record['ctime_ns']) == (ws['inode'], ws['mtime_ns'], ws['ctime_ns'])
        assert record['path'] == str(f.resolve()) and record['size'] == REAL_STAT(f).st_size


def test_hash_caches_miss_after_a_change(sealed):
    import structure_generation as SG
    f, ws = sealed
    first = SG._sha(f)
    assert SG._stat(f)['inode'] == ws['inode']
    f.write_bytes(b'xyz' * 1000)
    assert SG._stat(f)['inode'] != ws['inode']
    assert SG._sha(f) == SI.sha256(f) != first


WORKERS = '''
import multiprocessing as mp, os, sys
def ino(p):
    return os.stat(p).st_ino
if __name__ == "__main__":
    for method in ("fork", "spawn"):
        with mp.get_context(method).Pool(1) as pool:
            print(method, pool.apply(ino, (sys.argv[1],)))
'''


def test_fork_and_spawn_workers_keep_the_hook(tmp_path):
    f = tmp_path / 'a'
    f.write_text('x')
    seal, digest = make_seal(tmp_path, [f])
    script = tmp_path / 'workers.py'
    script.write_text(WORKERS)
    env = {k: v for k, v in os.environ.items() if not k.startswith('CODEX_SITE')}
    env.update(PYTHONPATH=str(HPG), CODEX_SITE_SEAL=str(seal), CODEX_SITE_SEAL_SHA256=digest)
    r = subprocess.run([sys.executable, str(script), str(f)], capture_output=True, text=True, env=env)
    ws = json.loads(seal.read_text())['files'][0]['workstation']['inode']
    assert r.returncode == 0, r.stderr
    assert r.stdout.split() == ['fork', str(ws), 'spawn', str(ws)]


def test_python_without_site_runs_without_the_hook(tmp_path):
    f = tmp_path / 'a'
    f.write_text('x')
    seal, digest = make_seal(tmp_path, [f])
    code = 'import os, sys; print(os.stat(sys.argv[1]).st_ino)'
    env = {'CODEX_SITE_SEAL': str(seal), 'CODEX_SITE_SEAL_SHA256': digest}
    on = run_python(code, env, str(f))
    off = run_python(code, env, str(f), flags=('-S',))
    assert int(on.stdout) == REAL_STAT(f).st_ino + 10 ** 6
    assert int(off.stdout) == REAL_STAT(f).st_ino


def test_journal_entry(monkeypatch, sealed):
    assert SI.journal_entry()['seal_sha256']
    monkeypatch.setattr(os, 'stat', REAL_STAT)
    monkeypatch.delenv('SLURM_JOB_ID', raising=False)
    assert SI.journal_entry() is None
    monkeypatch.setenv('SLURM_JOB_ID', '1')
    with pytest.raises(RuntimeError, match='must be active'):
        SI.journal_entry()


def seal_tool(*args):
    return subprocess.run([sys.executable, str(HPG / 'site_seal.py'), *map(str, args)], capture_output=True, text=True,
                          env={k: v for k, v in os.environ.items() if not k.startswith('CODEX_SITE')})


def test_capture_then_verify_and_seal(tmp_path):
    data = tmp_path / 'data'
    data.mkdir()
    files = [data / 'one.bin', data / 'with space.tif']
    for i, p in enumerate(files):
        p.write_bytes(bytes([i]) * 5000)
    (data / 'link').symlink_to(files[0])
    listing = tmp_path / 'list.txt'
    listing.write_text('\n'.join(map(str, [*files, data / 'link'])) + '\n')
    r = seal_tool('capture', '--list', listing, '--out', tmp_path / 'cap', '--workers', 2)
    assert r.returncode == 0, r.stderr
    digest = r.stdout.split()[-1]
    sums = (tmp_path / 'cap' / 'SHA256SUMS').read_text()
    assert str(files[1]) in sums and subprocess.run(['sha256sum', '-c', '--quiet', tmp_path / 'cap' / 'SHA256SUMS']).returncode == 0
    ids = tmp_path / 'cap' / 'workstation_identities.json'
    r = seal_tool('verify-and-seal', '--identities', ids, '--identities-sha256', digest, '--out', tmp_path / 'seal.json')
    assert r.returncode == 0, r.stderr
    seal = json.loads((tmp_path / 'seal.json').read_text())
    assert len(seal['files']) == 2 and seal['symlinks_checked'] == 1
    assert not os.access(tmp_path / 'seal.json', os.W_OK)
    table, _ = SI.load_table(tmp_path / 'seal.json', SI.sha256(tmp_path / 'seal.json'))
    assert set(table) == {REAL_STAT(p).st_ino for p in files}


def test_verify_refuses_a_changed_copy(tmp_path):
    f = tmp_path / 'a.bin'
    f.write_bytes(b'a' * 100)
    listing = tmp_path / 'list.txt'
    listing.write_text(f'{f}\n')
    r = seal_tool('capture', '--list', listing, '--out', tmp_path / 'cap')
    digest = r.stdout.split()[-1]
    f.write_bytes(b'b' * 100)
    r = seal_tool('verify-and-seal', '--identities', tmp_path / 'cap' / 'workstation_identities.json',
                  '--identities-sha256', digest, '--out', tmp_path / 'seal.json')
    assert r.returncode == 2 and not (tmp_path / 'seal.json').exists()
    failures = json.loads((tmp_path / 'seal.json.failures.json').read_text())['failures']
    assert 'sha256' in failures[0]['problem']
