#!/usr/bin/env python3
"""Refresh structures inside verified Tissue without rejecting edge-touching objects.

Original 15 retain their original classifiers. The five additions retain fitted Nerve/Capillary
thresholds and the original Lymphatics classifier. The QuPath stage never saves image data.
Run only after the candidate Tissue is installed and reviewed.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'typing'))
import codex_typing as CT
import structure_thresholds as ST

SCRIPT = HERE / 'qupath' / 'structures_in_tissue.groovy'
CLASSIFIERS = Path('/data/QuPath_Projects/QuPath_project_IsletExplorer/classifiers/pixel_classifiers')
QUPATH = '/home/smith6jt/QuPath/bin/QuPath'

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def threshold_inputs(run: Path, row: dict) -> dict | None:
    """Fingerprint the exact fitted values used for additions, without unrelated historical reports."""
    if row['cohort'] != 'addition':
        return None
    fitted = json.loads((run / 'structures' / 'thresholded.json').read_text())
    return {'reference_raw_threshold': fitted['reference_raw_threshold'], 'sigma_px': fitted['sigma_px']}

def one(run: Path, row: dict, force: bool = False) -> dict:
    tag = row['image_name']
    tissue = run / 'inputs' / 'tissue' / f'{tag}.geojson'
    project = run / 'qupath'
    entries = {e['imageName']: str(e['entryID']) for e in json.loads((project / 'project.qpproj').read_text())['images']}
    qpdata = project / 'data' / entries[tag] / 'data.qpdata'
    record = run / 'structures' / 'redraw_records' / f'{tag}.json'
    fingerprint = {'tissue_sha256': sha(tissue), 'qpdata_sha256': sha(qpdata), 'script_sha256': sha(SCRIPT),
                   'driver_sha256': sha(__file__), 'threshold_code_sha256': sha(ST.__file__),
                   'threshold_inputs': threshold_inputs(run, row),
                   'classifiers': {c: sha(CLASSIFIERS / f'{c}.json') for c in ('Nerve','Capillary','Lymphatics')}}
    if record.exists() and not force:
        old = json.loads(record.read_text())
        if old.get('inputs') == fingerprint and all(sha(run / 'structures' / c / f'{c}__{tag}.geojson') == r['sha256']
                                                   for c, r in old['classes'].items()):
            print(tag, 'structures already current', flush=True)
            return old
        raise ValueError(f'{tag}: stale structure checkpoint; use --force for candidate-only rebuild')
    before = fingerprint['qpdata_sha256']
    worker = run / 'work' / tag / ('structures_' + str(time.time_ns()))
    classifiers = worker / 'classifiers' / 'pixel_classifiers'
    classifiers.mkdir(parents=True)
    shutil.copy2(project / 'project.qpproj', worker / 'project.qpproj')
    (worker / 'data').symlink_to(project / 'data', target_is_directory=True)
    for c in fingerprint['classifiers']:
        shutil.copy2(CLASSIFIERS / f'{c}.json', classifiers / f'{c}.json')
    classes = ['Nerve','Capillary','Lymphatic'] if row['cohort'] == 'original' else ['Lymphatic']
    log = run / 'logs' / tag / 'structures_redraw.log'; log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [QUPATH, 'script', '-p', str(worker / 'project.qpproj'), '-i', tag,
           '-a', str(worker), '-a', ','.join(classes), '-a', str(tissue), str(SCRIPT)]
    with log.open('w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                env={**os.environ, 'JAVA_TOOL_OPTIONS':'-Xmx64g -XX:ActiveProcessorCount=12'})
    if result.returncode or any(not (worker / f'{c}.geojson').exists() for c in classes):
        raise RuntimeError(f'{tag}: structure classifier failed; {log}')
    if sha(qpdata) != before: raise RuntimeError(f'{tag}: read-only classifier changed image data')
    if row['cohort'] == 'addition':
        kept = run / 'structures' / 'superseded_additions' / 'Lymphatic'
        kept.mkdir(parents=True, exist_ok=True)
        shutil.copy2(worker / 'Lymphatic.geojson', kept / f'Lymphatic__{tag}.geojson')
        fitted = fingerprint['threshold_inputs']
        methods = ST.make_one(run, tag, fitted['reference_raw_threshold'], fitted['sigma_px'], retain_edge=True)
    else:
        methods = {}
        for c in classes:
            target = run / 'structures' / c / f'{c}__{tag}.geojson'; target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(worker / f'{c}.geojson', target)
            methods[c] = {'origin':'original classifier within current Tissue', 'edge_touching_retained':True}
    rec = {'image':tag, 'inputs':fingerprint, 'tissue_sha256':fingerprint['tissue_sha256'], 'classes':{}}
    for c in ('Nerve','Capillary','Lymphatic'):
        path = run / 'structures' / c / f'{c}__{tag}.geojson'
        rec['classes'][c] = {**methods[c], 'sha256':sha(path), 'objects':len(json.loads(path.read_text())['features']),
                             'edge_touching_retained':True}
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(json.dumps(rec, indent=2)+'\n')
    print(tag, {c:r['objects'] for c,r in rec['classes'].items()}, flush=True)
    return rec

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True); p.add_argument('--images', nargs='*')
    p.add_argument('--workers', type=int, default=2); p.add_argument('--force', action='store_true')
    a=p.parse_args(); run=a.run.resolve()
    rows=[r for r in CT.cohort_rows() if not a.images or r['image_name'] in a.images]
    with ThreadPoolExecutor(a.workers) as pool:
        records=list(pool.map(lambda r:one(run,r,a.force),rows))
    path=run/'structures'/'redrawn.json'
    old=json.loads(path.read_text()).get('images',{}) if path.exists() else {}
    old.update({r['image']:r for r in records})
    path.write_text(json.dumps({'images':old},indent=2)+'\n')
if __name__ == '__main__':main()
