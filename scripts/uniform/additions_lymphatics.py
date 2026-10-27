#!/usr/bin/env python3
"""Draw the additions' Lymphatic objects again inside their current Tissue (owner, 2026-10-27).

typing/structure_thresholds.py keeps the original `Lymphatics` classifier's objects for the five additions
(KEEP_CLASSIFIER), read from `$OUT/structures/superseded_additions/Lymphatic/`. Those were drawn in the review project
inside the defective Tissue outlines. This runs qupath/lymphatics_in_tissue.groovy read-only on each image's entry
of the uniform project (after tissue_import), with the original project's classifier file (the same file drew the
original 15's lymphatics), and replaces the kept objects; the earlier ones move to
`superseded_additions/Lymphatic_review_tissue/`. Run `structure_thresholds.py make` afterwards.

    python scripts/uniform/additions_lymphatics.py --out $OUT --images Aab_6090 ND_6232 T1D_6520 ND_6539 Aab_6558
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUPATH = '/home/smith6jt/QuPath/bin/QuPath'
CLASSIFIER = Path('/data/QuPath_Projects/QuPath_project_IsletExplorer/classifiers/pixel_classifiers/Lymphatics.json')
SCRIPT = HERE / 'qupath' / 'lymphatics_in_tissue.groovy'


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', required=True)
    a = ap.parse_args()
    out = a.out.resolve()
    project = out / 'qupath'
    entries = {e['imageName']: int(e['entryID']) for e in json.loads((project / 'project.qpproj').read_text())['images']}
    kept_dir = out / 'structures' / 'superseded_additions' / 'Lymphatic'
    earlier_dir = out / 'structures' / 'superseded_additions' / 'Lymphatic_review_tissue'
    records = {}
    for tag in a.images:
        tissue_rec = json.loads((out / 'provenance' / tag / 'tissue_import.json').read_text())
        if tissue_rec.get('status') != 'pass':
            raise SystemExit(f'{tag}: tissue_import has not passed')
        qpdata = project / 'data' / str(entries[tag]) / 'data.qpdata'
        before = sha256(qpdata)
        worker = out / 'work' / tag / f'{time.strftime("%Y%m%dT%H%M%S")}.lymphatics'
        (worker / 'classifiers' / 'pixel_classifiers').mkdir(parents=True)
        shutil.copy2(project / 'project.qpproj', worker / 'project.qpproj')
        (worker / 'data').symlink_to(project / 'data')
        shutil.copy2(CLASSIFIER, worker / 'classifiers' / 'pixel_classifiers' / 'Lymphatics.json')
        target = worker / f'Lymphatic__{tag}.geojson'
        log = out / 'logs' / tag / 'lymphatics.log'
        log.parent.mkdir(parents=True, exist_ok=True)
        cmd = [QUPATH, 'script', '-p', str(worker / 'project.qpproj'), '-i', tag, '-a', str(target), str(SCRIPT)]
        with open(log, 'a') as fh:
            fh.write(f'\n===== {time.strftime("%Y-%m-%dT%H:%M:%S%z")}\n{" ".join(cmd)}\n')
            fh.flush()
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                env={**os.environ, 'JAVA_TOOL_OPTIONS': '-Xmx96g'}).returncode
        if rc != 0 or not target.is_file():
            raise SystemExit(f'{tag}: QuPath failed ({log})')
        if sha256(qpdata) != before:
            raise SystemExit(f'{tag}: the image data changed during a read-only step')
        n = len(json.loads(target.read_text())['features'])
        earlier = kept_dir / f'Lymphatic__{tag}.geojson'
        if earlier.exists() and not (earlier_dir / earlier.name).exists():
            earlier_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(earlier), earlier_dir / earlier.name)
        kept_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target, kept_dir / f'Lymphatic__{tag}.geojson')
        records[tag] = {'objects': n, 'output': str(kept_dir / f'Lymphatic__{tag}.geojson'), 'sha256': sha256(target),
                        'tissue_import_attempt': tissue_rec.get('attempt_id'), 'log': str(log)}
        print(f'{tag}: {n} Lymphatic objects', flush=True)
    prov = out / 'structures' / 'additions_lymphatics.json'
    prev = json.loads(prov.read_text()) if prov.exists() else {}
    prev.update({'classifier': {'path': str(CLASSIFIER), 'sha256': sha256(CLASSIFIER)},
                 'script': {'path': str(SCRIPT), 'sha256': sha256(SCRIPT)}, 'time': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                 'images': {**prev.get('images', {}), **records}})
    prov.write_text(json.dumps(prev, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
