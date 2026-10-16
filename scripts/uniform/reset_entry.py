#!/usr/bin/env python3
"""Start images of the uniform project again from their first step: move everything run_uniform.py made for them
into <out>/superseded/<label>/<image>/, keeping relative paths. Nothing is deleted; source projects are not touched.

Use it when an image must go through the cell steps again after its islet steps ran (run_uniform.py refuses that,
because cell detection keeps existing annotations), or to keep a reviewed state aside before a full rerun.

    python scripts/uniform/reset_entry.py --out data/development/uniform-20donor-20261015 \\
        --label canary-v1-20261015 --images Aab_6505 ND_6232
"""
from __future__ import annotations

import argparse
import fcntl
import json
import shutil
from pathlib import Path

# Everything run_uniform.py writes per image ({TAG} = image name); directories are moved whole.
PER_IMAGE = (
    'provenance/{TAG}', 'logs/{TAG}', 'snapshots/{TAG}', 'work/{TAG}',
    'inputs/tissue/{TAG}.geojson', 'cells/{TAG}.parquet', 'planes/{TAG}.planes.tif', 'planes/{TAG}.planes.json',
    'cell_geojson/{TAG}.geojson', 'fit/{TAG}.json', 'fit/{TAG}.mask.npy', 'fit/{TAG}.objects.npy',
    'script1_outlines/{TAG}.geojson', 'islets/{TAG}.geojson', 'annotations/{TAG}.tsv',
    'final_cells/{TAG}.parquet', 'context/{TAG}.parquet', 'review/{TAG}.segmentation_comparison.json',
    'qupath/classifiers/pixel_classifiers/Islets_{TAG}.json',
    'qupath/classifiers/pixel_classifiers/IsletHormone_{TAG}.json',
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--label', required=True, help='name of the superseded set, e.g. canary-v1-20261015')
    ap.add_argument('--images', nargs='+', required=True)
    a = ap.parse_args()
    out = a.out.resolve()
    entries = {e['imageName']: int(e['entryID']) for e in json.loads((out / 'qupath' / 'project.qpproj').read_text())['images']}
    unknown = [t for t in a.images if t not in entries]
    if unknown:
        raise SystemExit(f'not in the uniform project: {unknown}')
    locks = []
    for tag in a.images:   # the same lock run_uniform.py holds while it works on an image
        (out / 'locks').mkdir(exist_ok=True)
        lock = open(out / 'locks' / f'{tag}.lock', 'w')
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit(f'{tag} is being run by run_uniform.py; not resetting it')
        locks.append(lock)
    for tag in a.images:
        dest = out / 'superseded' / a.label / tag
        if dest.exists():
            raise SystemExit(f'{dest} exists: choose another --label')
        moved = []
        entry_dir = out / 'qupath' / 'data' / str(entries[tag])
        for f in sorted(entry_dir.iterdir()) if entry_dir.is_dir() else []:   # the entry folder itself stays (empty)
            target = dest / 'qupath' / 'data' / str(entries[tag]) / f.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), target)
            moved.append(str(f.relative_to(out)))
        for pattern in PER_IMAGE:
            src = out / pattern.replace('{TAG}', tag)
            if src.exists():
                target = dest / pattern.replace('{TAG}', tag)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), target)
                moved.append(pattern.replace('{TAG}', tag))
        (dest / 'reset.json').write_text(json.dumps({'image': tag, 'entry_id': entries[tag], 'label': a.label,
                                                     'moved': moved}, indent=2) + '\n')
        print(f'{tag}: moved {len(moved)} item(s) to {dest}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
