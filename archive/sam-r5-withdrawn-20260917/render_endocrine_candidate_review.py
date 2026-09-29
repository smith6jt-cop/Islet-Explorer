#!/usr/bin/env python3
"""Render raw-channel checks at fixed image-wide ranges, with old/new boundaries."""
import argparse
import gzip
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts' / 'cohort_expansion'))
from prepare_endocrine_candidates import write_json
from convert_images import canonical_channel


def render(folder):
    spec = json.loads((folder/'prompts.json').read_text())
    pixels = tifffile.memmap(spec['source_image'], mode='r')
    with tifffile.TiffFile(spec['source_image']) as tf:
        labels = [canonical_channel(s) for s in tf.imagej_metadata['Labels']]
    with gzip.open(spec['source_geometry'], 'rt') as stream:
        old = {f['properties']['name']: shape(f['geometry']) for f in json.load(stream)['features']
               if f['properties'].get('classification', {}).get('name') == 'Islet'}
    new = {f['properties']['name']: shape(f['geometry'])
           for f in json.loads((folder/'review_candidates.geojson').read_text())['features']
           if f['geometry']['type'] in ['Polygon', 'MultiPolygon']}
    markers = ['INS', 'GCG', 'SST', 'CHGA']
    high = dict(zip(spec['rgb_channels'], spec['display_max']))
    low = dict(zip(spec['rgb_channels'], spec['display_min']))
    high['CHGA'] = float(np.quantile(pixels[labels.index('CHGA'), ::16, ::16], .998))
    low['CHGA'] = 0
    comparison = pd.read_csv(folder/'previous_islet_comparison.csv')
    positives = comparison.sort_values(['old_hormone_seed_cells', 'previous_islet'], ascending=[False, True]).head(6)
    if spec['case_id'] == '78':
        negatives = comparison.set_index('previous_islet').loc[
            ['Islet_130', 'Islet_206', 'Islet_1148', 'Islet_124', 'Islet_161', 'Islet_179']].reset_index()
    else:
        negatives = comparison.loc[comparison.old_hormone_seed_cells == 0].sort_values('old_area_um2').head(6)
    candidates = pd.read_csv(folder/'candidate_review.csv').sort_values(['hormone_seed_cells', 'candidate_id'])
    sampled = candidates.iloc[np.unique(np.linspace(0, len(candidates)-1, 9).round().astype(int))]
    selections = {
        'old_positive_controls': [(r.previous_islet, old[r.previous_islet], f'{r.old_hormone_seed_cells} old hormone seeds') for r in positives.itertuples()],
        'old_negative_controls': [(r.previous_islet, old[r.previous_islet], f'{r.old_hormone_seed_cells} old hormone seeds') for r in negatives.itertuples()],
        'new_candidates_size_range': [(r.candidate_id, new[r.candidate_id], f'{r.hormone_seed_cells} hormone seeds') for r in sampled.itertuples()],
    }
    old_tree, new_tree = shapely.STRtree(list(old.values())), shapely.STRtree(list(new.values()))
    records = []
    for title, examples in selections.items():
        if not examples: continue
        fig, axes = plt.subplots(len(examples), 5, figsize=(13, 2.5*len(examples)), squeeze=False)
        for row, (name, geom, note) in enumerate(examples):
            cx, cy = geom.centroid.x, geom.centroid.y
            half = max(110, int(max(geom.bounds[2]-geom.bounds[0], geom.bounds[3]-geom.bounds[1])/2)+40)
            x0, y0 = max(0, int(cx)-half), max(0, int(cy)-half)
            x1, y1 = min(pixels.shape[2], int(cx)+half), min(pixels.shape[1], int(cy)+half)
            box = shapely.box(x0, y0, x1, y1)
            rgb = np.stack([np.clip((pixels[labels.index(m), y0:y1, x0:x1].astype(float)-low[m])/(high[m]-low[m]), 0, 1)
                            for m in markers[:3]], axis=-1)
            for col, label in enumerate(['INS / GCG / SST']+markers):
                ax = axes[row, col]
                if col == 0: ax.imshow(rgb)
                else: ax.imshow(pixels[labels.index(label), y0:y1, x0:x1], cmap='gray', vmin=low[label], vmax=high[label])
                for tree, color in [(old_tree, '#ff9c24'), (new_tree, '#21e78b')]:
                    for idx in tree.query(box, predicate='intersects'):
                        for part in shapely.get_parts(tree.geometries[idx]):
                            if part.geom_type != 'Polygon': continue
                            xy = np.asarray(part.exterior.coords)
                            ax.plot(xy[:, 0]-x0, xy[:, 1]-y0, color=color, linewidth=.8)
                ax.set_xlim(0, x1-x0); ax.set_ylim(y1-y0, 0)
                ax.set_xticks([]); ax.set_yticks([])
                if row == 0: ax.set_title(label)
                if col == 0: ax.set_ylabel(name+'\n'+note, fontsize=8)
            records.append(dict(sheet=title, name=name, x_px=x0, y_px=y0, width_px=x1-x0, height_px=y1-y0))
        fig.suptitle(f'{spec["case_id"]}: orange = previous outlines; green = new unreviewed boundaries\n'
                     'Fixed per-image channel ranges. Point-only endocrine foci are available in QuPath.', fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, .95))
        for suffix in ['png', 'pdf']: fig.savefig(folder/f'{title}.{suffix}', dpi=140)
        plt.close(fig)
    write_json(folder/'raw_channel_review.json', dict(source_image=spec['source_image'], display_min=low,
        display_max=high, examples=records, raw_pixels=True, per_crop_contrast=False))
    print(f'{spec["case_id"]}: rendered {len(records)} raw-channel examples', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--case-ids', nargs='+', default=['78'])
    a = p.parse_args()
    for cid in a.case_ids: render(a.root/cid)


if __name__ == '__main__': main()
