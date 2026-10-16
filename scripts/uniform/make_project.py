#!/usr/bin/env python3
"""Create the uniform-pipeline QuPath project: one fresh entry per image in cohort.csv.

Each entry copies its source entry's server definition (image file, channel names, calibration), so it reads the
same pixels as the source project; it starts with no objects. The source projects are only read. The Tissue
outlines are imported afterwards by run_uniform.py (tissue_export on the source, tissue_import on this project).

    python scripts/uniform/make_project.py --out data/development/uniform-20donor-20261015
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import shutil
import sys
import time
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / 'scripts'))
from excluded_donors import assert_no_excluded  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def image_uri(builder: dict) -> str:
    while 'builder' in builder:
        builder = builder['builder']
    return builder['uri']


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True, help='fresh output root; the project goes to <out>/qupath')
    ap.add_argument('--cohort', type=Path, default=HERE / 'cohort.csv')
    args = ap.parse_args()
    out = args.out.resolve()
    project_dir = out / 'qupath'
    if project_dir.exists():
        raise SystemExit(f'{project_dir} exists: the uniform project is created once (use a fresh --out)')
    rows = list(csv.DictReader(args.cohort.open()))
    if len({r['image_name'] for r in rows}) != len(rows):
        raise SystemExit('duplicate image_name in the cohort')
    assert_no_excluded([r['case_id'] for r in rows], 'Uniform cohort')

    sources, entries, record = {}, [], []
    for new_id, row in enumerate(rows, start=1):
        src_dir = Path(row['source_project'])
        if not src_dir.is_absolute():
            src_dir = REPO / src_dir
        src_file = src_dir / 'project.qpproj'
        if src_file not in sources:
            sources[src_file] = json.loads(src_file.read_text())
        match = [e for e in sources[src_file]['images'] if int(e['entryID']) == int(row['source_entry_id'])]
        if len(match) != 1:
            raise SystemExit(f"{row['image_name']}: entry {row['source_entry_id']} not found once in {src_file}")
        src = match[0]
        uri = image_uri(src['serverBuilder'])
        image_path = Path(uri.replace('file://', '', 1) if uri.startswith('file://') else uri.replace('file:', '', 1))
        image_path = Path(str(image_path).replace('%20', ' '))
        if not image_path.is_file():
            raise SystemExit(f"{row['image_name']}: image file missing: {image_path}")
        entries.append({
            'serverBuilder': copy.deepcopy(src['serverBuilder']),
            'entryID': new_id,
            'randomizedName': str(uuid.uuid4()),
            'imageName': row['image_name'],
            'metadata': {'Case ID': row['case_id'], 'Donor status': row['donor_status'], 'Cohort': row['cohort'],
                         'Source project': str(src_dir), 'Source entry': str(row['source_entry_id']),
                         'Pipeline': 'CODEX uniform: InstanSeg 0.1.1 on 9 channels; Senior v3 islets (CHGA)'},
            'tags': [],
        })
        st = image_path.stat()
        record.append({'entry_id': new_id, **row, 'source_project_file': str(src_file),
                       'source_project_sha256': sha256(src_file), 'source_image_name': src['imageName'],
                       'image_path': str(image_path), 'image_bytes': st.st_size, 'image_mtime_ns': st.st_mtime_ns})

    now = int(time.time() * 1000)
    project = {'version': '0.7.0', 'createTimestamp': now, 'modifyTimestamp': now,
               'uri': (project_dir / 'project.qpproj').as_uri(), 'lastID': len(entries), 'images': entries}
    project_dir.mkdir(parents=True)
    for e in entries:
        (project_dir / 'data' / str(e['entryID'])).mkdir(parents=True)
    (project_dir / 'classifiers' / 'pixel_classifiers').mkdir(parents=True)
    (project_dir / 'scripts' / 'lib').mkdir(parents=True)
    shutil.copy2(HERE / 'qupath' / 'lib' / 'common.groovy', project_dir / 'scripts' / 'lib' / 'common.groovy')
    (project_dir / 'project.qpproj').write_text(json.dumps(project, indent=2) + '\n')
    manifest = {'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'project': str(project_dir / 'project.qpproj'),
                'cohort_file': str(args.cohort.resolve()), 'cohort_sha256': sha256(args.cohort),
                'common_groovy_sha256': sha256(project_dir / 'scripts' / 'lib' / 'common.groovy'),
                'entries': record}
    (out / 'project_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'{len(entries)} entries -> {project_dir / "project.qpproj"}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
