#!/usr/bin/env python3
"""Run the CODEX uniform chain on images of the uniform project (owner directives 2026-10-15).

Every image gets the same steps: its existing Tissue outline is copied from its source project, cells and nuclei
are segmented with InstanSeg fluorescence_nuclei_and_cells 0.1.1 on 9 channels, then Islet-Explorer-Senior's v3
islet method is applied (per-image classifier with CHGA in CD99's place, script 1, step 1b on the cell outlines).

QuPath steps run on a private copy of the project file (QuPath never writes the project file itself) with Senior's
footer appended (qupath/footer_template.groovy). The footer checks the step's postcondition, saves the image data
only when it passed (save policy footer-saves; step 0 saves itself after its own checks), and writes the provenance
JSON that is the step's success record. A step whose prerequisites have no passing record is refused.

    python scripts/uniform/run_uniform.py --out data/development/uniform-20donor-20261015 \\
        --images Aab_6505 ND_6232 [--steps tissue_export tissue_import ...] [--force]
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
QUPATH = '/home/smith6jt/QuPath/bin/QuPath'
FOOTER = HERE / 'qupath' / 'footer_template.groovy'
CFG_PLACEHOLDER = '/*__V3_CFG__*/ null'
GPU_ENV = {'PYTORCH_VERSION': '2.7.1', 'PYTORCH_FLAVOR': 'cu128',
           'PYTORCH_LIBRARY_PATH': '/home/smith6jt/.djl.ai/pytorch/2.7.1-cu128-linux-x86_64'}
GPU_JAVA = '-Dai.djl.offline=true -Dai.djl.pytorch.num_interop_threads=1 -Dai.djl.pytorch.num_threads=1'


@dataclass(frozen=True)
class Step:
    id: str
    kind: str                       # source (QuPath on the source project) | qupath | fit
    script: str = ''                # relative to scripts/uniform
    args: tuple = ()
    save: str = 'never'             # footer-saves | self-saves | never
    heap_gb: int = 64
    threads: int = 8
    timeout_min: int = 60
    postcondition: str = 'v3_result'
    after: tuple = ()
    gpu: bool = False
    post: str = ''                  # cells_parquet | context_parquet | ''
    overrides: tuple = field(default_factory=tuple)


ISLET_CLASSES = 'Islet|EndocrineCluster|EndocrineSingle'
STEPS = [
    Step('tissue_export', 'source', 'qupath/export_tissue.groovy', ('{OUT}/inputs/tissue/{TAG}.geojson',), timeout_min=30),
    Step('tissue_import', 'qupath', 'qupath/import_tissue.groovy', ('{OUT}/inputs/tissue/{TAG}.geojson',), 'footer-saves',
         timeout_min=30, after=('tissue_export',)),
    Step('preflight', 'qupath', 'qupath/preflight.groovy', (), 'footer-saves', timeout_min=15, postcondition='preflight',
         after=('tissue_import',)),
    Step('0', 'qupath', 'qupath/0_cell_detection.groovy', ('{CFG}',), 'self-saves', 96, 12, 180,
         'regions_detected(classes=Tissue)', after=('preflight',), gpu=True),
    Step('2', 'qupath', 'qupath/2_cell_measurements.groovy', (), 'footer-saves', 160, 24, 240,
         'cells_measured(measurement=Cell: DAPI: Mean)', after=('0',), overrides=(('N_THREADS', '24'),)),
    Step('C_cells', 'qupath', 'qupath/exports.groovy', ('cells', '{TMP}/{TAG}.cells.tsv'), heap_gb=128, after=('2',),
         post='cells_parquet'),
    Step('C_planes', 'qupath', 'qupath/export_classifier_planes.groovy', ('{OUT}/planes/{TAG}.planes.tif', '2.0'),
         heap_gb=96, after=('tissue_import',)),
    Step('C_geojson', 'qupath', 'qupath/export_cells_geojson.groovy', ('{OUT}/cell_geojson/{TAG}.geojson',), heap_gb=96,
         after=('0',)),
    Step('fit', 'fit', 'islet_classifier_fit.py', after=('C_cells', 'C_planes')),
    Step('1', 'qupath', 'qupath/1_islet_objects.groovy', (), 'footer-saves', 96, 16, 240,
         f'islets_present(classes={ISLET_CLASSES})', after=('fit',)),
    Step('1_outlines', 'qupath', 'qupath/exports.groovy', ('islets', '{OUT}/script1_outlines/{TAG}.geojson'), heap_gb=96,
         after=('1',)),
    Step('1b', 'qupath', 'qupath/1b_islet_outlines_to_cells.groovy', (), 'footer-saves', 96, 16, 120,
         after=('1', '0', '1_outlines')),
    Step('4', 'qupath', 'qupath/4_islet_naming.groovy', (), 'footer-saves', 64, 16, 120,
         'islet_bands(classes=Islet;bands=none)', after=('1b',)),
    Step('3', 'qupath', 'qupath/3_signed_distance.groovy', (), 'footer-saves', 64, 16, 120,
         f'signed_distance(classes={ISLET_CLASSES}|Hormone Fragment|Artifact|pLN|Tissue)', after=('0', '1b')),
    Step('8', 'qupath', 'qupath/8_hierarchy_parenting.groovy', (), 'footer-saves', 64, 16, 120,
         'parenting(forbid=ROOT|NONE)', after=('0', '4')),
    Step('J_cells', 'qupath', 'qupath/exports.groovy', ('cells', '{TMP}/{TAG}.final_cells.tsv'), heap_gb=128,
         after=('3', '8'), post='cells_parquet'),
    Step('J_context', 'qupath', 'qupath/exports.groovy', ('context', '{TMP}/{TAG}.context.tsv'), heap_gb=96,
         after=('3', '8'), post='context_parquet'),
    Step('J_annotations', 'qupath', 'qupath/exports.groovy', ('annotations', '{OUT}/annotations/{TAG}.tsv'), heap_gb=96,
         after=('8',)),
    Step('J_islets', 'qupath', 'qupath/exports.groovy', ('islets', '{OUT}/islets/{TAG}.geojson'), heap_gb=96,
         after=('4',)),
]
BY_ID = {s.id: s for s in STEPS}
POST_TARGET = {('C_cells', 'cells_parquet'): 'cells/{TAG}.parquet', ('J_cells', 'cells_parquet'): 'final_cells/{TAG}.parquet',
               ('J_context', 'context_parquet'): 'context/{TAG}.parquet'}


def dependents(step_id: str) -> list[str]:
    """Every step that rests, directly or through others, on `step_id` (in STEPS order)."""
    found, frontier = set(), {step_id}
    while frontier:
        frontier = {s.id for s in STEPS if set(s.after) & frontier} - found
        found |= frontier
    return [s.id for s in STEPS if s.id in found]


# The islet chain starts from the image as these steps left it; a copy is kept (snapshots/<tag>/before_islets.qpdata)
# and restored before every run of step 1, so re-running the islets is exact rather than incremental.
CELL_STEPS = ('tissue_import', 'preflight', '0', '2')
ISLET_CHAIN = tuple(dependents('fit'))
# Steps that need every cohort image's cells (the cohort scale) run after all images have finished the others.
COHORT_PHASE = ('fit', *ISLET_CHAIN)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def stat_of(path: Path) -> dict:
    try:
        st = os.stat(path)
    except FileNotFoundError:
        return {'exists': False, 'path': str(path)}
    return {'exists': True, 'path': str(path), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}


def parse_postcondition(text: str) -> dict:
    """Senior v3ctl.parse_postcondition: name(key=v1|v2;key2=v)."""
    m = re.fullmatch(r'\s*([A-Za-z_]\w*)\s*(?:\((.*)\))?\s*', text)
    if not m:
        raise ValueError(f'bad postcondition {text!r}')
    params = {}
    for part in (m.group(2) or '').split(';'):
        if part.strip():
            key, _, val = part.partition('=')
            params[key.strip()] = [x.strip() for x in val.split('|')]
    return {'name': m.group(1), 'params': params}


def groovy_literal(v) -> str:
    """Senior v3ctl.groovy_literal: a Groovy literal for a JSON-like value."""
    if v is None:
        return 'null'
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, int):
        return f'{v}L' if abs(v) > 2 ** 31 - 1 else str(v)
    if isinstance(v, float):
        return repr(v) + 'd'
    if isinstance(v, (str, Path)):
        out = []
        for ch in str(v):
            if ch == '\\':
                out.append('\\\\')
            elif ch == "'":
                out.append("\\'")
            elif ch == '\n':
                out.append('\\n')
            elif ch == '\r':
                out.append('\\r')
            elif ch == '\t':
                out.append('\\t')
            elif ord(ch) < 32 or ord(ch) > 126:
                units = ch.encode('utf-16-be')
                out += ['\\u%04X' % int.from_bytes(units[i:i + 2], 'big') for i in range(0, len(units), 2)]
            else:
                out.append(ch)
        return "'" + ''.join(out) + "'"
    if isinstance(v, (list, tuple)):
        return '[' + ', '.join(groovy_literal(x) for x in v) + ']'
    if isinstance(v, dict):
        if not v:
            return '[:]'
        return '[' + ', '.join(f'{groovy_literal(str(k))}: {groovy_literal(x)}' for k, x in v.items()) + ']'
    raise TypeError(f'no Groovy literal for {type(v)}')


def apply_overrides(body: str, overrides: tuple) -> tuple[str, dict]:
    """Senior steps.tsv overrides: rewrite the declaration '<type> NAME = ...' (must match exactly once)."""
    applied = {}
    for name, value in overrides:
        pattern = re.compile(rf'^(\s*)((?:int|double|boolean|String|long|float)\s+{re.escape(name)})\s*=\s*[^\n]*$', re.M)
        hits = pattern.findall(body)
        if len(hits) != 1:
            raise SystemExit(f'override {name}: expected one declaration, found {len(hits)}')
        body = pattern.sub(lambda m: f'{m.group(1)}{m.group(2)} = {value}', body, count=1)
        applied[name] = value
    return body, applied


def git_head() -> str:
    try:
        return subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    except OSError:
        return ''


class Runner:
    def __init__(self, out: Path, rows: list[dict], cohort_path: Path):
        self.out = out
        self.project_dir = out / 'qupath'
        self.project = self.project_dir / 'project.qpproj'
        if not self.project.is_file():
            raise SystemExit(f'{self.project} missing: run make_project.py first')
        entries = json.loads(self.project.read_text())['images']
        self.entry_of = {e['imageName']: int(e['entryID']) for e in entries}
        self.rows = {r['image_name']: r for r in rows}
        self.cohort_path = cohort_path
        self.run_id = time.strftime('%Y%m%dT%H%M%S') + f'-{os.getpid()}'
        self.tmp = out / 'tmp'

    def record(self, tag: str, step: str) -> Path:
        return self.out / 'provenance' / tag / f'{step}.json'

    def passed(self, tag: str, step: str) -> bool:
        f = self.record(tag, step)
        if not f.is_file():
            return False
        try:
            rec = json.loads(f.read_text())
        except json.JSONDecodeError:
            return False
        if rec.get('status') != 'pass':
            return False
        if step == 'fit':   # stale if made with another cohort scale, if any image's cells changed since the scale,
            import cohort_scale  # noqa: PLC0415   if the scale was built under another rule, or by other fit code
            made_by = rec.get('script') or {}
            if (made_by.get('sha256'), made_by.get('library_sha256')) != (sha256(HERE / 'islet_classifier_fit.py'),
                                                                          sha256(HERE / 'islet_cell_test.py')):
                return False
            scale = self.cohort_scale_file()
            used = (rec.get('cohort_scale') or {}).get('sha256')
            if not (used and scale.is_file() and sha256(scale) == used):
                return False
            current = json.loads(scale.read_text())
            return cohort_scale.follows_rule(current) and not cohort_scale.holds_current(current, self.out, list(self.rows))
        return True

    def cohort_scale_file(self) -> Path:
        return self.out / 'fit' / 'cohort_scale.json'

    def qpdata_of(self, tag: str) -> Path:
        return self.project_dir / 'data' / str(self.entry_of[tag]) / 'data.qpdata'

    def snapshot_of(self, tag: str) -> Path:
        return self.out / 'snapshots' / tag / 'before_islets.qpdata'

    def retire(self, path: Path, note: str) -> None:
        """Move a record or snapshot aside to history/ (never deleted)."""
        if path.exists():
            hist = path.parent / 'history'
            hist.mkdir(parents=True, exist_ok=True)
            path.rename(hist / f'{path.stem}.{time.time_ns()}.{note}{path.suffix}')

    def invalidate(self, tag: str, step_id: str) -> None:
        """Running a step makes every step that rests on it stale; a cell step also retires the pre-islet snapshot."""
        stale = [s for s in dependents(step_id) if self.record(tag, s).exists()]
        for s in stale:
            self.retire(self.record(tag, s), f'stale-after-{step_id}')
        if stale:
            print(f'{tag}: {step_id}: records of {stale} moved to history (they rest on this step)', flush=True)
        if step_id in CELL_STEPS and self.snapshot_of(tag).exists():
            self.retire(self.snapshot_of(tag), f'stale-after-{step_id}')
            self.retire(self.snapshot_of(tag).with_suffix('.json'), f'stale-after-{step_id}')
            print(f'{tag}: {step_id}: pre-islet snapshot retired', flush=True)

    def start_islets(self, tag: str) -> bool:
        """Before step 1: restore the image as the cell steps left it, or take that copy on the first islet run."""
        snap, qp = self.snapshot_of(tag), self.qpdata_of(tag)
        meta = snap.with_suffix('.json')
        if snap.is_file():
            m = json.loads(meta.read_text())
            if sha256(snap) != m.get('sha256'):
                print(f'{tag}: 1: REFUSED, the pre-islet snapshot does not match its record ({meta})', flush=True)
                return False
            part = qp.with_name(qp.name + '.restore')
            shutil.copy2(snap, part)
            os.replace(part, qp)
            print(f'{tag}: 1: restored the pre-islet snapshot taken {m.get("taken")}', flush=True)
            return True
        prov = self.out / 'provenance' / tag
        touched = [s for s in ISLET_CHAIN if self.record(tag, s).exists() or any((prov / 'history').glob(f'{s}.*.json'))]
        if touched:
            print(f'{tag}: 1: REFUSED, islet steps already ran on this entry ({touched}) and no pre-islet snapshot exists; '
                  f'reset the entry (README: "Starting an image again")', flush=True)
            return False
        snap.parent.mkdir(parents=True, exist_ok=True)
        part = snap.with_name(snap.name + '.part')
        shutil.copy2(qp, part)
        digest = sha256(part)
        part.rename(snap)
        cell_records = {s: json.loads(self.record(tag, s).read_text()).get('attempt_id') for s in CELL_STEPS
                        if self.record(tag, s).is_file()}
        meta.write_text(json.dumps({'source': str(qp), 'sha256': digest, 'size': snap.stat().st_size,
                                    'taken': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'cell_step_attempts': cell_records},
                                   indent=2) + '\n')
        print(f'{tag}: 1: pre-islet snapshot taken ({snap})', flush=True)
        return True

    def ensure_cohort_scale(self) -> Path | None:
        """The cohort scale holding every cohort image's current cells under the current rule; (re)built when every
        image's cells exist."""
        import cohort_scale  # noqa: PLC0415  (beside this script)
        target = self.cohort_scale_file()
        tags = list(self.rows)
        if target.is_file():
            current = json.loads(target.read_text())
            if cohort_scale.follows_rule(current) and not cohort_scale.holds_current(current, self.out, tags):
                return target
        not_ready = [t for t in tags if not self.passed(t, 'C_cells')]
        if not_ready:
            print(f'fit: REFUSED, the cohort scale needs every image\'s cells; C_cells has not passed for {not_ready}', flush=True)
            return None
        if target.is_file():
            self.retire(target, 'rebuilt')
        log = self.out / 'logs' / '_cohort' / 'cohort_scale.log'
        log.parent.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(HERE / 'cohort_scale.py'), '--out', str(self.out), '--cohort', str(self.cohort_path)]
        with open(log, 'a') as fh:
            fh.write(f'\n===== {time.strftime("%Y-%m-%dT%H:%M:%S%z")}\n{" ".join(cmd)}\n')
            fh.flush()
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL).returncode
        ok = rc == 0 and target.is_file()
        print(f'cohort scale: {"built" if ok else "FAILED"} ({log})', flush=True)
        return target if ok else None

    def subst(self, text: str, tag: str, cfg_path: Path | None = None) -> str:
        return (text.replace('{OUT}', str(self.out)).replace('{TMP}', str(self.tmp)).replace('{TAG}', tag)
                .replace('{CFG}', str(cfg_path) if cfg_path else ''))

    def run(self, tag: str, step: Step, force: bool) -> bool:
        if self.passed(tag, step.id) and not force:
            print(f'{tag}: {step.id}: already passed (use --force to rerun)', flush=True)
            return True
        missing = [a for a in step.after if not self.passed(tag, a)]
        if missing:
            print(f'{tag}: {step.id}: REFUSED, prerequisites without a passing record: {missing}', flush=True)
            return False
        if step.id == '1' and not self.start_islets(tag):
            return False
        self.invalidate(tag, step.id)
        rec = self.record(tag, step.id)
        if rec.exists():
            hist = rec.parent / 'history'
            hist.mkdir(parents=True, exist_ok=True)
            rec.rename(hist / f'{step.id}.{time.time_ns()}.json')
        t0 = time.time()
        if step.id == 'tissue_export' and self.rows[tag].get('tissue_outline'):
            ok = self.copy_tissue_outline(tag, step)
        else:
            ok = self.run_fit(tag, step) if step.kind == 'fit' else self.run_qupath(tag, step)
        if ok and step.post:
            ok = self.post(tag, step)
        print(f'{tag}: {step.id}: {"PASS" if ok else "FAIL"} ({time.time() - t0:.0f}s)', flush=True)
        return ok

    # ---- QuPath steps ------------------------------------------------------------------------------------------
    def run_qupath(self, tag: str, step: Step) -> bool:
        row = self.rows[tag]
        if step.kind == 'source':
            src = Path(row['source_project'])
            project_dir = src if src.is_absolute() else REPO / src
            entry_id = int(row['source_entry_id'])
            src_entries = json.loads((project_dir / 'project.qpproj').read_text())['images']
            image_name = [e['imageName'] for e in src_entries if int(e['entryID']) == entry_id][0]
        else:
            project_dir, entry_id, image_name = self.project_dir, self.entry_of[tag], tag
        qpdata = project_dir / 'data' / str(entry_id) / 'data.qpdata'
        entry_dir = qpdata.parent
        dir_before = {f.name: stat_of(f) for f in entry_dir.iterdir()} if step.kind == 'source' else None
        attempt = f'{self.run_id}.{step.id}'
        worker = self.out / 'work' / tag / attempt
        worker.mkdir(parents=True)
        shutil.copy2(project_dir / 'project.qpproj', worker / 'project.qpproj')
        for d in ('data', 'classifiers', 'scripts', 'resources'):
            if (project_dir / d).exists():
                (worker / d).symlink_to(project_dir / d)
        src_script = HERE / step.script
        body, applied = apply_overrides(src_script.read_text(encoding='utf-8'), step.overrides)
        cfg_path = worker / 'step_config.json'
        args = [self.subst(a, tag, cfg_path) for a in step.args]
        for a in args:
            if a.startswith(str(self.out)):
                Path(a).parent.mkdir(parents=True, exist_ok=True)
        self.tmp.mkdir(parents=True, exist_ok=True)
        script_path = worker / 'step.groovy'
        before = stat_of(qpdata)
        cfg = {
            'run_id': self.run_id, 'attempt_id': attempt, 'seq': 1, 'step': step.id, 'passes': step.kind,
            'image': image_name, 'tag': tag, 'entry_id': entry_id, 'project': str(project_dir / 'project.qpproj'),
            'site': socket.gethostname(), 'worker_project': str(worker / 'project.qpproj'), 'out': str(self.out),
            'provenance': str(self.record(tag, step.id)), 'save_policy': step.save,
            'postcondition': parse_postcondition(step.postcondition), 'lane': 'gpu' if step.gpu else 'cpu', 'slot': None,
            'heap_gb': step.heap_gb, 'threads': step.threads, 'host': socket.gethostname(),
            'gpu': 'gpu0' if step.gpu else None, 'source_script': str(src_script), 'source_sha256': sha256(src_script),
            'script_selection': 'codex-uniform', 'repo_git': {'head': git_head()},
            'body_sha256': hashlib.sha256(body.encode('utf-8')).hexdigest(), 'instantiated_script': str(script_path),
            'overrides': applied, 'args': args, 'footer_template_sha256': sha256(FOOTER),
            'qpdata_before': before, 'driver_start': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'step_config': str(cfg_path),
            'deadline_epoch_ms': None, 'save_margin_min': 10.0,
        }
        cfg_path.write_text(json.dumps(cfg, indent=2) + '\n')
        footer = FOOTER.read_text(encoding='utf-8')
        if footer.count(CFG_PLACEHOLDER) != 1:
            raise SystemExit(f'{FOOTER}: expected exactly one {CFG_PLACEHOLDER!r}')
        script_path.write_text(body.rstrip('\n') + '\n' + footer.replace(CFG_PLACEHOLDER, groovy_literal(cfg)),
                               encoding='utf-8')
        env = dict(os.environ)
        java = f'-Xmx{step.heap_gb}g -Djava.util.concurrent.ForkJoinPool.common.parallelism={step.threads}'
        if step.gpu:
            env.update(GPU_ENV)
            java += ' ' + GPU_JAVA
        env['JAVA_TOOL_OPTIONS'] = java
        cmd = [QUPATH, 'script', '-p', str(worker / 'project.qpproj'), '-i', image_name]
        for a in args:
            cmd += ['-a', a]
        cmd.append(str(script_path))
        log = self.out / 'logs' / tag / f'{step.id}.log'
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, 'a') as fh:
            fh.write(f'\n===== {time.strftime("%Y-%m-%dT%H:%M:%S%z")} {attempt}\n{" ".join(cmd)}\n')
            fh.flush()
            try:
                rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, stdin=subprocess.DEVNULL,
                                    timeout=step.timeout_min * 60).returncode
            except subprocess.TimeoutExpired:
                fh.write(f'TIMEOUT after {step.timeout_min} min\n')
                rc = -1
        if step.kind == 'source':
            dir_after = {f.name: stat_of(f) for f in entry_dir.iterdir()}
            if dir_after != dir_before:
                changed = sorted(k for k in set(dir_before) | set(dir_after) if dir_before.get(k) != dir_after.get(k))
                print(f'{tag}: {step.id}: SOURCE entry files changed during a read-only step: {entry_dir} {changed}')
                return False
        return rc == 0 and self.passed(tag, step.id)

    # ---- Python steps -----------------------------------------------------------------------------------------
    def copy_tissue_outline(self, tag: str, step: Step) -> bool:
        """tissue_export for an image whose cohort row names a redrawn outline (redraw_tissue.py, owner 2026-10-27):
        the file is copied in place of the source project's Tissue, and the record keeps both checksums."""
        src = REPO / self.rows[tag]['tissue_outline']
        if not src.is_file():
            print(f'{tag}: {step.id}: REFUSED, redrawn outline missing: {src}', flush=True)
            return False
        feats = json.loads(src.read_text())['features']
        # A named annotation changes every cell's `Parent` from "Annotation (Tissue)", which the pipeline keys tissue
        # cells on (islet_cell_test.TISSUE_PARENT): the first redraw import of 2026-10-27 left all cells ineligible.
        if len(feats) != 1 or (feats[0].get('properties') or {}).get('name') or \
                ((feats[0].get('properties') or {}).get('classification') or {}).get('name') != 'Tissue':
            print(f'{tag}: {step.id}: REFUSED, {src} must hold one unnamed Tissue annotation', flush=True)
            return False
        dst = Path(self.subst(step.args[0], tag))
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        ok = sha256(src) == sha256(dst)
        rec = {'schema': 'uniform.tissue_outline/1', 'status': 'pass' if ok else 'fail', 'step': step.id, 'tag': tag,
               'passes': 'redrawn', 'source': {'path': str(src), 'sha256': sha256(src)},
               'output': {'path': str(dst), 'sha256': sha256(dst)},
               'script': {'path': str(HERE / 'redraw_tissue.py'), 'sha256': sha256(HERE / 'redraw_tissue.py')},
               'repo_git': {'head': git_head()}, 'time': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
        self.record(tag, step.id).parent.mkdir(parents=True, exist_ok=True)
        self.record(tag, step.id).write_text(json.dumps(rec, indent=2) + '\n')
        return ok

    def run_fit(self, tag: str, step: Step) -> bool:
        scale = self.ensure_cohort_scale()
        if scale is None:
            return False
        planes = self.out / 'planes' / f'{tag}.planes.tif'
        cells = self.out / 'cells' / f'{tag}.parquet'
        classifier = self.project_dir / 'classifiers' / 'pixel_classifiers' / f'Islets_{tag}.json'
        report = self.out / 'fit' / f'{tag}.json'
        report.parent.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(HERE / step.script), '--tag', tag, '--planes', str(planes), '--cells', str(cells),
               '--cohort-scale', str(scale), '--classifier', str(classifier), '--report', str(report),
               '--mask', str(report.with_suffix('.mask.npy')), '--objects', str(report.with_suffix('.objects.npy'))]
        log = self.out / 'logs' / tag / f'{step.id}.log'
        log.parent.mkdir(parents=True, exist_ok=True)
        started = time.strftime('%Y-%m-%dT%H:%M:%S%z')
        with open(log, 'a') as fh:
            fh.write(f'\n===== {started}\n{" ".join(cmd)}\n')
            fh.flush()
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL).returncode
        marker = classifier.with_name('IsletHormone_' + classifier.name[len('Islets_'):])
        ok = rc == 0 and report.is_file() and classifier.is_file() and marker.is_file()
        rec = {'status': 'pass' if ok else 'fail', 'step': step.id, 'tag': tag, 'command': cmd, 'exit_code': rc,
               'started': started, 'ended': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
               'script': {'path': str(HERE / step.script), 'sha256': sha256(HERE / step.script),
                          'library_sha256': sha256(HERE / 'islet_cell_test.py')},
               'inputs': {str(p): sha256(p) for p in (planes, cells) if p.is_file()},
               'cohort_scale': {'file': str(scale), 'sha256': sha256(scale)},
               'outputs': {str(p): sha256(p) for p in (classifier, marker, report) if p.is_file()}}
        self.record(tag, step.id).parent.mkdir(parents=True, exist_ok=True)
        self.record(tag, step.id).write_text(json.dumps(rec, indent=2) + '\n')
        return ok

    def post(self, tag: str, step: Step) -> bool:
        """TSV export -> Parquet, checked against the export's own detection count; the TSV is then removed."""
        import pyarrow.csv as pcsv
        import pyarrow.parquet as pq
        tsv = Path(self.subst(step.args[1], tag))
        target = self.out / POST_TARGET[(step.id, step.post)].replace('{TAG}', tag)
        target.parent.mkdir(parents=True, exist_ok=True)
        table = pcsv.read_csv(tsv, parse_options=pcsv.ParseOptions(delimiter='\t'),
                              convert_options=pcsv.ConvertOptions(strings_can_be_null=True))
        rec_path = self.record(tag, step.id)
        prov = json.loads(rec_path.read_text())
        expected = (prov.get('counts') or {}).get('detections')
        ids = table.column('Object ID').to_pandas()
        problems = []
        if expected is not None and table.num_rows != expected:
            problems.append(f'{table.num_rows} rows, {expected} detections')
        if ids.isna().any() or ids.duplicated().any():
            problems.append('missing or duplicate Object ID')
        if step.post == 'cells_parquet' and 'Cell: DAPI: Mean' not in table.column_names:
            problems.append("no 'Cell: DAPI: Mean' column")
        if problems:
            print(f'{tag}: {step.id}: post {step.post} FAILED: {problems}')
            prov['status'] = 'fail'
            prov['post'] = {'hook': step.post, 'problems': problems}
            rec_path.write_text(json.dumps(prov, indent=2) + '\n')
            return False
        part = target.with_suffix('.parquet.part')
        pq.write_table(table, part, compression='zstd')
        part.rename(target)
        tsv.unlink()
        prov['post'] = {'hook': step.post, 'output': str(target), 'rows': table.num_rows, 'columns': table.num_columns,
                        'sha256': sha256(target)}
        rec_path.write_text(json.dumps(prov, indent=2) + '\n')
        return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--cohort', type=Path, default=HERE / 'cohort.csv')
    ap.add_argument('--images', nargs='+', required=True, help='image names (e.g. Aab_6505) or "all"')
    ap.add_argument('--steps', nargs='+', default=[s.id for s in STEPS], choices=[s.id for s in STEPS])
    ap.add_argument('--force', action='store_true', help='rerun steps that already passed')
    args = ap.parse_args()
    out = args.out.resolve()
    rows = list(csv.DictReader(args.cohort.open()))
    sys.path.insert(0, str(REPO / 'scripts'))
    from excluded_donors import assert_no_excluded
    assert_no_excluded([r['case_id'] for r in rows], 'Uniform cohort')
    names = [r['image_name'] for r in rows] if args.images == ['all'] else args.images
    unknown = [n for n in names if n not in {r['image_name'] for r in rows}]
    if unknown:
        raise SystemExit(f'not in the cohort: {unknown}')
    order = [s for s in STEPS if s.id in args.steps]
    # One lock per image (as Senior's run_v3): each image writes only its own entry and classifier files, so
    # different images may run at the same time; the same image never runs twice at once.
    (out / 'locks').mkdir(exist_ok=True)
    locks = []
    for tag in names:
        lock = open(out / 'locks' / f'{tag}.lock', 'w')
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit(f'{tag} is being run by another run_uniform.py ({out / "locks" / f"{tag}.lock"})')
        locks.append(lock)
    runner = Runner(out, rows, args.cohort.resolve())
    failed = {}
    # Two phases: every image's own steps first, then the fit and the islet chain, which need the cohort scale and
    # therefore every cohort image's cells.
    for phase in ([s for s in order if s.id not in COHORT_PHASE], [s for s in order if s.id in COHORT_PHASE]):
        for tag in names:
            if tag in failed:
                continue
            for step in phase:
                if not runner.run(tag, step, args.force):
                    failed[tag] = step.id
                    break
    print('FAILED: ' + ', '.join(f'{t}:{s}' for t, s in failed.items()) if failed else 'all requested steps passed', flush=True)
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
