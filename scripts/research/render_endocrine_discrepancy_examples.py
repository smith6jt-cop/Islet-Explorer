#!/usr/bin/env python3
"""Render traced real cells with original object (cyan) and final cell (white).

Panels contain median-filtered image pixels, not a REDSEA-corrected image.
Numerical measurements and stage attribution live in the companion traces.
"""
import argparse
import json
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from shapely.geometry import shape
import typing_accuracy_tiles_20261029 as renderer

EXAMPLES = {
    'ND_0112': ('35550667-3158-4724-a948-34b513ec8455', 'A: positive mean clamped to zero by REDSEA'),
    'Aab_6558': ('1d783132-731e-4697-8967-9e32f358c3e8', 'B: localized classifier signal; whole-cell mean below level'),
    'T1D_6551': ('adcc76bf-ec08-4ff5-8b8a-693c71d1c870', 'C: tiny overlap expands to a whole-cell object'),
}


def features(path):
    raw = json.loads(path.read_text())
    return raw['features'] if isinstance(raw, dict) else raw


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--classifier-run', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        p.error('Choose a fresh output directory')
    a.out.mkdir(parents=True)
    renderer.shown_markers = lambda *_: (['INS', 'GCG', 'SST', 'CHGA'], [])
    done, records = [], []
    for tag, (oid, description) in EXAMPLES.items():
        d = pd.read_parquet(a.run / 'typing/identity' / f'{tag}.parquet')
        selected = d[d.object_id.eq(oid)].assign(image=tag, group='discrepancy', pair=description)
        if len(selected) != 1:
            raise ValueError(f'{tag}/{oid}: expected exactly one cell')
        row = renderer.render_image(tag, selected, a.run)[0]
        cell = next(shape(f['geometry']) for f in features(a.run / 'cell_geojson' / f'{tag}.geojson') if f['id'] == oid)
        original = [f for f in features(a.classifier_run / 'script1_outlines' / f'{tag}.geojson')
                    if f.get('properties', {}).get('classification', {}).get('name') == 'Islet']
        original = max(original, key=lambda f: cell.intersection(shape(f['geometry'])).area)
        original_geometry = shape(original['geometry'])
        q = pd.read_parquet(a.run / 'cells' / f'{tag}.parquet',
                            columns=['Object ID', 'Centroid X µm', 'Centroid Y µm']).set_index('Object ID').loc[oid]
        entry = renderer.CT.manifest(a.run)[tag]
        img = renderer.RS.Channels(entry['image_path'])
        H, W = img.shape[1:]
        img.close()
        half = round(renderer.CROP_UM / renderer.CT.PIXEL_UM / 2)
        cx, cy = q['Centroid X µm'] / .3774, q['Centroid Y µm'] / .3774
        origin = np.array([int(max(0, min(W - 2 * half, cx - half))),
                           int(max(0, min(H - 2 * half, cy - half)))])
        for components in renderer.polygon_components(original['geometry']):
            for ring in components:
                points = np.round((ring - origin) * renderer.ZOOM).astype(np.int32)
                for tile in row['tiles']:
                    cv2.polylines(tile, [points], True, (0, 255, 255), 1)
        records.append({'image': tag, 'cell_id': oid, 'example': description,
                        'original_object_id': original['id'],
                        'original_object_cell_overlap_fraction': cell.intersection(original_geometry).area / cell.area})
        done.append(row)
    renderer.sheets(done, a.out)
    (a.out / 'examples.json').write_text(json.dumps({'legend': __doc__, 'examples': records}, indent=2) + '\n')


if __name__ == '__main__':
    main()
