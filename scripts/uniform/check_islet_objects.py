#!/usr/bin/env python3
"""Check that every endocrine object is made of whole cell segmentations (owner, 2026-10-16).

The owner's rules for the objects step 1b leaves (Islet, EndocrineCluster, EndocrineSingle):
  - an object's outline follows the cell segmentations it picks up: no cell is cut by an object's edge;
  - no cell outline holds more than one object, so a single cell is never counted twice or as a cluster;
  - no endocrine object is smaller than a cell: every object holds at least one whole cell.

Read from the exports, independently of the QuPath script: `islets/<image>.geojson` (the objects) and
`cell_geojson/<image>.geojson` (the cells; a cell is a detection with a nucleus, as in step 1b). A cell is IN an object
when at least MEMBER_SHARE of its area lies inside; it is CUT when a share between EDGE_SHARE and 1 - EDGE_SHARE lies
inside (slivers below EDGE_SHARE are the cell outlines' and the objects' own pixel staircases).

    python scripts/uniform/check_islet_objects.py --out data/development/uniform-20donor-20261015 --images all
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely import STRtree
from shapely.geometry import shape

HERE = Path(__file__).resolve().parent
ENDOCRINE = ('Islet', 'EndocrineCluster', 'EndocrineSingle')
FRAGMENT = 'Hormone Fragment'
MEMBER_SHARE = 0.5
EDGE_SHARE = 0.05


def features(path: Path) -> list:
    raw = json.loads(path.read_text())
    return raw['features'] if isinstance(raw, dict) else raw


def valid(g):
    return g if g.is_valid else shapely.make_valid(g)


def check(out: Path, tag: str) -> dict:
    objs, classes = [], []
    n_fragments = 0
    for f in features(out / 'islets' / f'{tag}.geojson'):
        cls = ((f.get('properties') or {}).get('classification') or {}).get('name')
        if cls == FRAGMENT:
            n_fragments += 1
        if cls in ENDOCRINE or cls == FRAGMENT:
            objs.append(valid(shape(f['geometry'])))
            classes.append(cls)
    cells = [valid(shape(f['geometry'])) for f in features(out / 'cell_geojson' / f'{tag}.geojson')
             if f.get('nucleusGeometry') and f.get('geometry')]
    cells = np.asarray(cells, dtype=object)
    cell_area = shapely.area(cells)
    o_idx, c_idx = STRtree(cells).query(np.asarray(objs, dtype=object), predicate='intersects')
    share = shapely.area(shapely.intersection(np.asarray(objs, dtype=object)[o_idx], cells[c_idx])) / cell_area[c_idx]
    keep = share > EDGE_SHARE
    o_idx, c_idx, share = o_idx[keep], c_idx[keep], share[keep]
    members = np.bincount(o_idx[share >= MEMBER_SHARE], minlength=len(objs))
    holders = np.bincount(c_idx, minlength=len(cells))
    classes = np.asarray(classes)
    cut = (share < 1 - EDGE_SHARE)
    return {
        'image': tag, 'objects': len(objs), 'hormone_fragments': n_fragments,
        'objects_without_a_whole_cell': int((members == 0).sum()),
        'cells_in_two_or_more_objects': int((holders >= 2).sum()),
        'cells_cut_by_an_object_edge': int(np.unique(c_idx[cut]).size),
        'singles_not_one_cell': int(((classes == 'EndocrineSingle') & (members != 1)).sum()),
        'clusters_or_islets_with_one_cell_or_none': int((np.isin(classes, ['EndocrineCluster', 'Islet']) & (members < 2)).sum()),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', required=True, help='image names or "all"')
    ap.add_argument('--cohort', type=Path, default=HERE / 'cohort.csv')
    ap.add_argument('--csv', type=Path, default=None, help='also write the table here')
    a = ap.parse_args()
    tags = [r['image_name'] for r in csv.DictReader(a.cohort.open())] if a.images == ['all'] else a.images
    rows = []
    for tag in tags:
        if not (a.out / 'islets' / f'{tag}.geojson').is_file():
            print(f'{tag}: no islet objects exported', flush=True)
            continue
        rows.append(check(a.out, tag))
        print(json.dumps(rows[-1]), flush=True)
    if a.csv:
        with a.csv.open('w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    bad = [r['image'] for r in rows if any(r[k] for k in r if k not in ('image', 'objects'))]
    print(f'{len(rows) - len(bad)} of {len(rows)} images pass' + (f'; failing: {bad}' if bad else ''))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
