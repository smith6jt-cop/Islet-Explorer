#!/usr/bin/env python3
"""Audit final SAM geometry against raw hormone pixels and render a pilot baseline.

This writes diagnostics only. It neither replaces native objects nor accepts
the provisional pixel cutoffs as a biological definition of an islet.
"""
import argparse
from datetime import datetime, timezone
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
import tifffile

from endocrine_signal_support import (hormone_support, mask_support_metrics,
                                      signal_footprint, keep_seed_near_components)
from prepare_endocrine_candidates import digest, write_json


def raw_crop(pixels, spec, geom, pad_um=12):
    pad = int(np.ceil(pad_um/spec['pixel_size_um']))
    a, b, c, d = geom.bounds
    x0, y0 = max(0, int(np.floor(a))-pad), max(0, int(np.floor(b))-pad)
    x1, y1 = min(pixels.shape[2], int(np.ceil(c))+pad+1), min(pixels.shape[1], int(np.ceil(d))+pad+1)
    raw = np.stack([pixels[spec['source_channel_names'].index(m), y0:y1, x0:x1]
                    for m in spec['rgb_channels']], axis=-1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    mask = shapely.intersects_xy(geom, xx, yy)
    return raw, mask, (x0, y0, x1, y1)


def audit_donor(root, destination, cid, fraction, tolerance):
    folder = root/cid
    out = destination/cid
    out.mkdir(exist_ok=False)
    spec = json.loads((folder/'prompts.json').read_text())
    image = Path(spec['source_image']).resolve(strict=True)
    previous = json.loads((root/'source-image-verification.json').read_text())
    source = next(s for s in previous['images'] if s['case_id'] == cid)
    stat = image.stat()
    if str(image) != source['path'] or stat.st_size != source['bytes'] or stat.st_mtime_ns != source['mtime_ns']:
        raise ValueError('Source image fingerprint differs')
    overlay = folder/'review_candidates.geojson'
    expected = json.loads((folder/'candidate_summary.json').read_text())['overlay_sha256']
    if digest(overlay) != expected:
        raise ValueError('Preserve changed candidate geometry for separate review')
    features = {f['properties']['name']: f for f in json.loads(overlay.read_text())['features']
                if f['geometry']['type'] in ('Polygon', 'MultiPolygon')}
    rows = []
    pixels = tifffile.memmap(image, mode='r')
    thresholds = [spec['operational_thresholds'][m] for m in spec['rgb_channels']]
    for name, feature in features.items():
        raw, mask, _ = raw_crop(pixels, spec, shape(feature['geometry']))
        metrics = mask_support_metrics(mask, hormone_support(raw, thresholds, spec['pixel_size_um']),
                                       spec['pixel_size_um'], tolerance)
        measures = feature['properties']['measurements']
        rows.append(dict(case_id=cid, candidate_id=name,
                         classification=feature['properties']['classification']['name'],
                         area_um2=measures['Candidate: area um2'],
                         hormone_seed_cells=measures['Candidate: hormone-positive seed cells'],
                         total_cells=measures['Candidate: total cells in new boundary'],
                         **metrics,
                         low_signal_proximity=metrics['area_within_signal_distance_fraction'] < fraction))
    table = pd.DataFrame(rows)
    table.to_csv(out/'mask_support.csv', index=False)
    # Include the actual reported failure, worst green masks, and strong/large controls.
    selected = ['New_78_0242'] if cid == '78' else []
    green = table[table.classification.str.startswith('Islet candidate')]
    selected += green.sort_values(['area_within_signal_distance_fraction', 'candidate_id']).head(2).candidate_id.tolist()
    selected += table.sort_values(['area_um2', 'candidate_id'], ascending=[False, True]).head(2).candidate_id.tolist()
    selected += table.sort_values(['area_within_signal_distance_fraction', 'candidate_id'], ascending=[False, True]).head(1).candidate_id.tolist()
    selected = list(dict.fromkeys(selected))
    seed_table = pd.read_csv(folder/'seed_cells.csv')
    seed_xy = seed_table[['X_centroid', 'Y_centroid']].to_numpy()/spec['pixel_size_um']
    display = json.loads((folder/'raw_channel_review.json').read_text())
    fig, axes = plt.subplots(len(selected), 4, figsize=(13, 3*len(selected)), squeeze=False)
    examples = []
    for index, name in enumerate(selected):
        feature = features[name]
        geom = shape(feature['geometry'])
        raw, mask, bounds = raw_crop(pixels, spec, geom)
        x0, y0, x1, y1 = bounds
        support = hormone_support(raw, thresholds, spec['pixel_size_um'])
        footprint = signal_footprint(support, spec['pixel_size_um'])
        nearby = seed_xy[((seed_xy >= [x0, y0]) & (seed_xy < [x1, y1])).all(axis=1)]-[x0, y0]
        footprint = keep_seed_near_components(footprint, nearby, spec['pixel_size_um'])
        rgb = np.clip((raw-np.asarray(spec['display_min']))/
                      (np.asarray(spec['display_max'])-spec['display_min']), 0, 1)
        chga = pixels[spec['source_channel_names'].index('CHGA'), y0:y1, x0:x1]
        row = table.set_index('candidate_id').loc[name]
        for col, title in enumerate(['Raw INS/GCG/SST', 'Current SAM boundary', 'Signal footprint pilot', 'CHGA context']):
            ax = axes[index, col]
            if col == 3:
                ax.imshow(chga, cmap='gray', vmin=display['display_min']['CHGA'], vmax=display['display_max']['CHGA'])
            else:
                ax.imshow(rgb)
            if col == 1:
                for part in shapely.get_parts(geom):
                    for ring in [part.exterior, *part.interiors]:
                        xy = np.asarray(ring.coords)-[x0, y0]
                        ax.plot(xy[:, 0], xy[:, 1], color='yellow', linewidth=.7)
            elif col == 2 and footprint.any() and not footprint.all():
                ax.contour(footprint.astype(float), levels=[.5], colors=['cyan'], linewidths=.7)
            ax.set_xlim(0, x1-x0); ax.set_ylim(y1-y0, 0)
            ax.set_xticks([]); ax.set_yticks([])
            if index == 0: ax.set_title(title)
            if col == 0:
                ax.set_ylabel(f'{name}\n{row.area_within_signal_distance_fraction:.0%} within {tolerance:g} um of signal\n{row.hormone_seed_cells}/{row.total_cells} seed cells', fontsize=8)
        examples.append(dict(candidate_id=name, source_pixel_bounds=bounds,
                             old_mask_pixels=int(mask.sum()), footprint_pixels=int(footprint.sum()),
                             footprint_does_not_define_whole_islet=True))
    fig.suptitle(f'Donor {cid}: SAM failure audit and conservative signal-footprint pilot\n'
                 'Fixed display ranges. Pixel cutoffs are provisional; neither outline is accepted for analysis.', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .96))
    for ext in ['png', 'pdf']: fig.savefig(out/f'comparison.{ext}', dpi=160)
    plt.close(fig)
    summary = dict(case_id=cid, boundaries=len(table),
                   flagged_low_signal_proximity=int(table.low_signal_proximity.sum()),
                   green_boundaries=len(green),
                   green_flagged_low_signal_proximity=int(green.low_signal_proximity.sum()),
                   examples=examples, input_overlay_sha256=digest(overlay),
                   source_image_fingerprint=source,
                   diagnostic_cutoffs=spec['operational_thresholds'],
                   cutoffs_validated_for_pixel_segmentation=False)
    write_json(out/'summary.json', summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('examples','source_image_fingerprint','diagnostic_cutoffs')}), flush=True)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate-root', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--case-ids', nargs='+', default=['78','6090','6232','6520','6539','6558'])
    p.add_argument('--minimum-proximity-fraction', type=float, default=.8)
    p.add_argument('--distance-um', type=float, default=3)
    a = p.parse_args()
    if not 0 < a.minimum_proximity_fraction <= 1 or a.distance_um <= 0:
        raise ValueError('Invalid diagnostic criterion')
    root = a.candidate_root.resolve(strict=True)
    out = a.output_root.resolve()
    out.mkdir(parents=True, exist_ok=False)
    summaries = [audit_donor(root,out,cid,a.minimum_proximity_fraction,a.distance_um) for cid in a.case_ids]
    record = dict(status='diagnostic_pilot_only', created_at_utc=datetime.now(timezone.utc).isoformat(),
                  candidate_root=str(root), donors=summaries,
                  minimum_proximity_fraction=a.minimum_proximity_fraction, distance_um=a.distance_um,
                  smoothing_sigma_um=.6, footprint_closing_radius_um=2, footprint_max_distance_um=3,
                  footprint_max_filled_hole_area_um2=80, footprint_seed_proximity_um=4,
                  native_objects_modified=False, analytical_acceptance=False,
                  limitations=['Cell-mean thresholds reused as provisional pixel cutoffs; calibration remains unvalidated.',
                               'Signal footprints can fragment whole islets and cannot identify hormone-negative endocrine structures.',
                               'A proximity flag is a review priority, not a biological rejection or acceptance.'],
                  code=[dict(path=str(path),sha256=digest(path)) for path in [Path(__file__),Path(__file__).with_name('endocrine_signal_support.py')]])
    write_json(out/'audit_summary.json', record)
    content = ''.join(f'<h2>Donor {html.escape(s["case_id"])}</h2>'
                      f'<p>{s["flagged_low_signal_proximity"]}/{s["boundaries"]} masks flagged by the provisional proximity check; '
                      f'{s["green_flagged_low_signal_proximity"]}/{s["green_boundaries"]} previously green masks flagged.</p>'
                      f'<a href="{s["case_id"]}/comparison.pdf"><img src="{s["case_id"]}/comparison.png"></a>' for s in summaries)
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>SAM boundary quality audit</title>'
        '<style>body{font:16px system-ui;max-width:1400px;margin:2rem auto;padding:1rem}img{max-width:100%}</style>'
        '<h1>SAM boundary failure audit</h1><p>R5 is withdrawn as a usable boundary result. Yellow is its SAM boundary; '
        'cyan is a conservative raw-signal footprint pilot. Pixel thresholds are provisional. No replacement islet boundaries '
        'have been accepted or written to QuPath.</p>'+content+'</html>')


if __name__ == '__main__': main()
