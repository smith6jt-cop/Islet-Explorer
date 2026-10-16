#!/usr/bin/env python3
"""The cohort's median p99.9 of every identity marker, over all images of the uniform cohort (2026-10-15).

islet_classifier_fit.py reads a cell's highest identity marker with each marker divided by the image's own p99.9
(islet_cell_test.common_scale). Insulin or somatostatin (islet_cell_test.FLOORED_HORMONES, recorded here as
`floor_markers`) with a p99.9 in the image below the cohort median written here is divided by the median instead:
an image can lack their cells almost entirely (beta cells in T1D, delta cells in many images), and its own p99.9 is
then background. Every cohort image's cells Parquet must exist; the file records each one's size, mtime and SHA-256.
islet_classifier_fit.py refuses a scale that does not hold its image's current cells or was built under another
rule (quantile, floor markers), and run_uniform.py rebuilds it then.

    python scripts/uniform/cohort_scale.py --out data/development/uniform-20donor-20261015
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import islet_cell_test as ict  # noqa: E402

SCALE_FILE = 'fit/cohort_scale.json'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def cells_path(out: Path, tag: str) -> Path:
    return out / 'cells' / f'{tag}.parquet'


def file_id(path: Path) -> dict:
    st = path.stat()
    return {'size': st.st_size, 'mtime_ns': st.st_mtime_ns}


def rule() -> dict:
    """What a scale is built under; a scale built under anything else is stale."""
    return {'quantile': ict.SCALE_QUANTILE, 'floor_markers': list(ict.FLOORED_HORMONES)}


def follows_rule(scale: dict) -> bool:
    return all(scale.get(k) == v for k, v in rule().items())


def holds_current(scale: dict, out: Path, tags: list[str]) -> list[str]:
    """Images whose cells are missing from `scale` or changed since it was built (empty: the scale is current)."""
    images = scale.get('images', {})
    stale = [t for t in tags if t not in images or not cells_path(out, t).is_file()
             or {k: images[t].get(k) for k in ('size', 'mtime_ns')} != file_id(cells_path(out, t))]
    return stale + sorted(set(images) - set(tags))


def build(out: Path, tags: list[str]) -> dict:
    import islet_classifier_fit as fit  # noqa: PLC0415  (cv2 and friends only when a scale is built)
    missing = [t for t in tags if not cells_path(out, t).is_file()]
    if missing:
        raise SystemExit(f'the cohort scale needs every image\'s cells; missing: {missing}')
    images = {}
    for tag in tags:
        path = cells_path(out, tag)
        _, values, _, _, eligible, *_ = fit.load_cells(path, None, None)
        images[tag] = {'cells_file': str(path), **file_id(path), 'sha256': sha256(path),
                       'eligible_cells': int(eligible.sum()), 'tops': ict.image_tops(values, eligible)}
        print(f'{tag}: ' + ', '.join(f'{h} {images[tag]["tops"].get(h, float("nan")):.0f}' for h in ict.HORMONES), flush=True)
    markers = sorted(set.intersection(*(set(v['tops']) for v in images.values())))
    return {
        'schema': 1, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), **rule(),
        'median': {m: float(np.nanmedian([images[t]['tops'][m] for t in tags])) for m in markers},
        'script': {'path': str(Path(__file__).resolve()), 'sha256': sha256(Path(__file__).resolve()),
                   'library_sha256': sha256(HERE / 'islet_cell_test.py')},
        'images': images,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--cohort', type=Path, default=HERE / 'cohort.csv')
    a = ap.parse_args()
    out = a.out.resolve()
    tags = [r['image_name'] for r in csv.DictReader(a.cohort.open())]
    scale = build(out, tags)
    target = out / SCALE_FILE
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix('.json.part')
    part.write_text(json.dumps(scale, indent=2) + '\n')
    part.rename(target)
    print(json.dumps({'file': str(target), 'images': len(tags),
                      'median': {h: round(scale['median'][h], 1) for h in ict.HORMONES}}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
