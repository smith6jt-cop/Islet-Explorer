#!/usr/bin/env python3
"""Systematic disease-blinded image panels and geometric checks for rebuilt islets.

Rendering and geometric checks do not imply biological acceptance. Actual visual
review decisions are recorded in an evidence-bound, separate JSON file. Source
TIFFs are read-only memory maps or three-channel reads; original pixels never change.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('MPLCONFIGDIR', '/tmp/islet-boundary-audit-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import distance_transform_edt
import shapely
from shapely.geometry import shape
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cohort_expansion'))
from verify_size_adaptive_output import verify
from common import sha256, write_json

MARKERS = ['INS', 'GCG', 'SST']
SEED = 20260920


def contour(ax, geom, offset, color='cyan'):
    for part in shapely.get_parts(geom):
        for ring in [part.exterior, *part.interiors]:
            xy = np.asarray(ring.coords) - offset
            ax.plot(xy[:, 0], xy[:, 1], color=color, linewidth=.7)


def select_examples(table):
    """Fixed categories independent of disease and image appearance."""
    chosen, seen = [], set()
    def add(frame, reason, mode='first'):
        frame = frame.loc[~frame.candidate_id.isin(seen)]
        if not len(frame):
            return
        r = frame.iloc[len(frame) // 2 if mode == 'median' else 0].to_dict()
        r['selection_reason'] = reason
        chosen.append(r)
        seen.add(r['candidate_id'])
    by_area = table.sort_values(['area_um2', 'candidate_id'])
    add(by_area.loc[by_area.area_um2 < 2000], 'small: median section area', 'median')
    add(by_area.loc[by_area.area_um2.between(2000, 10000, inclusive='left')], 'medium: median section area', 'median')
    add(by_area.loc[by_area.area_um2 >= 10000].sort_values('area_um2', ascending=False), 'large: largest section area')
    add(table.sort_values(['foreground_fraction', 'candidate_id']), 'lowest classifier foreground fraction')
    add(table.sort_values(['max_gap_to_foreground_um', 'candidate_id'], ascending=[False, True]), 'largest internal foreground gap')
    order = table.candidate_id.map(lambda s: hashlib.sha256(f'{SEED}\0{s}'.encode()).hexdigest())
    random = table.assign(selection_hash=order).sort_values('selection_hash')
    add(random, 'deterministic random')
    while len(chosen) < min(6, len(table)):
        add(random, 'deterministic random: fills absent size stratum')
    return pd.DataFrame(chosen)


def support_metrics(features, foreground, pixel_um, downsample):
    records = []
    for f in features:
        geom = shape(f['geometry'])
        name = f['properties']['name']
        a, b, c, d = geom.bounds
        x0, x1 = max(0, int(np.floor(a / downsample)) - 2), min(foreground.shape[1], int(np.ceil(c / downsample)) + 3)
        y0, y1 = max(0, int(np.floor(b / downsample)) - 2), min(foreground.shape[0], int(np.ceil(d / downsample)) + 3)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        region = shapely.intersects_xy(geom, (xx + .5) * downsample, (yy + .5) * downsample)
        crop = foreground[y0:y1, x0:x1]
        distances = distance_transform_edt(~crop) * downsample * pixel_um
        measures = f['properties'].get('measurements', {})
        fraction = measures.get('Adaptive boundary: foreground_fraction')
        records.append({'candidate_id': name, 'area_um2': geom.area * pixel_um ** 2,
            'foreground_fraction': float(fraction) if fraction is not None else (float(crop[region].mean()) if region.any() else np.nan),
            'max_gap_to_foreground_um': float(distances[region].max()) if region.any() else np.nan,
            'gap_gt_3um_fraction': float((distances[region] > 3).mean()) if region.any() else np.nan,
            'gap_gt_10um_fraction': float((distances[region] > 10).mean()) if region.any() else np.nan,
            'sampled_region_pixels': int(region.sum()), 'valid_geometry': bool(geom.is_valid)})
    return pd.DataFrame(records)


def audit_one(boundaries_root, output_root, cid, adaptive_subdir='adaptive'):
    boundaries_root, output_root = Path(boundaries_root), Path(output_root)
    source_file = boundaries_root / 'sources' / cid / 'source.json'
    adaptive = boundaries_root / adaptive_subdir / cid
    source = json.loads(source_file.read_text())
    prediction_file = adaptive / 'prediction_manifest.json'
    prediction = json.loads(prediction_file.read_text())
    geometry_file = adaptive / 'review_candidates.geojson'
    fingerprint = {'source_specification_sha256': sha256(source_file), 'geometry_sha256': sha256(geometry_file),
                   'prediction_manifest_sha256': sha256(prediction_file)}
    out = output_root / cid
    out.mkdir(parents=True, exist_ok=True)
    existing = out / 'render_manifest.json'
    if existing.exists():
        record = json.loads(existing.read_text())
        if record['input_fingerprint'] != fingerprint:
            raise ValueError(f'Inputs changed since boundary audit for {cid}; use a fresh output root')
        if not (out / 'review_panel.png').is_file() or sha256(out / 'review_panel.png') != record['panel_sha256']:
            raise ValueError(f'Boundary review panel changed or missing: {cid}')
        return record
    if prediction['markers'] != MARKERS or int(prediction['operating_point']['minimum_votes']) != 160:
        raise ValueError('Boundary classifier differs from fixed three-marker 160-vote specification')
    image_file = Path(source['source_image'])
    stat = image_file.stat()
    if stat.st_size != source['source_bytes'] or stat.st_mtime_ns != source['source_mtime_ns']:
        raise ValueError(f'Changed source image fingerprint for {cid}')
    geometry_check = verify(adaptive)
    features = json.loads(geometry_file.read_text())['features']
    by_name = {f['properties']['name']: shape(f['geometry']) for f in features}
    foreground_file = boundaries_root / 'predictions' / cid / 'pixel_classes.png'
    foreground = np.asarray(Image.open(foreground_file)) > 0
    px = float(source['pixel_size_um'])
    downsample = float(prediction['prediction_downsample'])
    table = support_metrics(features, foreground, px, downsample)
    table.to_csv(out / 'all_boundary_metrics.csv', index=False)
    selected = select_examples(table)
    selected.to_csv(out / 'selected_examples.csv', index=False)
    channel_indices = [source['source_channel_names'].index(m) for m in MARKERS]
    try:
        image = tifffile.memmap(image_file, mode='r')
        image_reader = 'read_only_memory_map'
    except ValueError:
        # Some source TIFFs are tiled/compressed. Decode the three requested
        # channels only; never materialize all 31 channels or change the TIFF.
        image = tifffile.imread(image_file, key=channel_indices, maxworkers=1)
        channel_indices = list(range(3))
        image_reader = 'three_channel_read_only_decode'
    if list(image.shape[1:]) != source['image_shape_cyx'][1:]:
        raise ValueError('Unexpected source image axis order')
    # Fixed per-donor display from classifier's tissue background normalization;
    # no per-crop rescaling. Values are for visualization, never positivity calls.
    low = np.asarray(prediction['normalization']['offset'], dtype=float)
    high = low + np.asarray(prediction['normalization']['scale'], dtype=float) * 12
    fig, axes = plt.subplots(3, 4, figsize=(16, 12), squeeze=False)
    rendered = []
    all_geometries = list(by_name.values())
    tree = shapely.STRtree(all_geometries)
    for i, record in enumerate(selected.to_dict('records')):
        geom = by_name[record['candidate_id']]
        a, b, c, d = geom.bounds
        pad = max(12 / px, .08 * max(c - a, d - b))
        x0, x1 = max(0, int(np.floor(a - pad))), min(image.shape[2], int(np.ceil(c + pad)))
        y0, y1 = max(0, int(np.floor(b - pad))), min(image.shape[1], int(np.ceil(d + pad)))
        raw = np.stack([image[ch, y0:y1, x0:x1] for ch in channel_indices], axis=-1).astype(float)
        rgb = np.clip((raw - low) / np.maximum(high - low, 1), 0, 1)
        row, start_col = divmod(i, 2)
        start_col *= 2
        for variant in range(2):
            ax = axes[row, start_col + variant]
            ax.imshow(rgb, interpolation='nearest')
            if variant:
                for j in tree.query(shapely.box(x0, y0, x1, y1), predicate='intersects'):
                    contour(ax, all_geometries[j], np.array([x0, y0]), '#808080')
                contour(ax, geom, np.array([x0, y0]), 'cyan')
            ax.set_xlim(0, x1 - x0)
            ax.set_ylim(y1 - y0, 0)
            ax.set_xticks([])
            ax.set_yticks([])
            label = record['candidate_id'].rsplit('_', 1)[-1]
            ax.set_title(f'{label}: {record["selection_reason"]}\n{record["area_um2"]:.0f} µm²; foreground {record["foreground_fraction"]:.0%}; gap {record["max_gap_to_foreground_um"]:.1f} µm' if not variant else 'Same raw field + selected outline', fontsize=8)
            # Physical scale bar retained independently of crop dimensions.
            bar_um = 10 if max(c - a, d - b) * px < 100 else 50
            bar_px = bar_um / px
            if bar_px < (x1 - x0) * .6:
                ax.plot([.06 * (x1-x0), .06 * (x1-x0) + bar_px], [.94 * (y1-y0)] * 2, color='white', lw=2)
                ax.text(.06 * (x1-x0), .90 * (y1-y0), f'{bar_um} µm', color='white', fontsize=7)
        rendered.append({'candidate_id': record['candidate_id'], 'selection_reason': record['selection_reason'],
                         'source_pixel_bounds': [x0, y0, x1, y1], 'review_status': 'unreviewed'})
    for i in range(len(selected) * 2, 12):
        axes.flat[i].set_axis_off()
    fig.suptitle(f'Image reference {cid} — systematic regenerated-boundary review\nINS red; GCG green; SST blue. Fixed within-image display. No disease category shown.\nCyan: selected boundary; gray: nearby boundaries. Geometric foreground fraction is not biological identity.', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .94))
    for ext in ['png', 'pdf']:
        fig.savefig(out / f'review_panel.{ext}', dpi=120)
    plt.close(fig)
    record = {'case_id': cid, 'status': 'rendered_awaiting_visual_review', 'input_fingerprint': fingerprint,
        'objects': len(features), 'source_image_modified': False, 'source_image_hash_recomputed': False, 'image_reader': image_reader,
        'geometry_check': geometry_check, 'selected_examples': rendered,
        'panel_sha256': sha256(out / 'review_panel.png'),
        'display': {'markers_rgb': MARKERS, 'minimum': low.tolist(), 'maximum': high.tolist(), 'per_crop_rescaling': False},
        'sampling': 'Six deterministic examples: small median, medium median, largest, lowest foreground fraction, largest foreground gap, random; missing strata filled by random.',
        'objects_foreground_below_80pct': int((table.foreground_fraction < .8).sum()),
        'objects_internal_gap_gt10um': int((table.max_gap_to_foreground_um > 10).sum()),
        'visual_review_completed': False, 'owner_approval_claimed': False}
    write_json(existing, record)
    print(json.dumps({'case_id': cid, 'objects': len(features), 'panel': str(out / 'review_panel.png')}), flush=True)
    return record


def render_context(boundaries_root, output_root, requests, adaptive_subdir='adaptive'):
    """Render unresolved examples with independent nuclear/endocrine context."""
    boundaries_root, output_root = Path(boundaries_root), Path(output_root)
    result = []
    for cid, names in requests.items():
        source = json.loads((boundaries_root / 'sources' / cid / 'source.json').read_text())
        prediction = json.loads((boundaries_root / adaptive_subdir / cid / 'prediction_manifest.json').read_text())
        features = json.loads((boundaries_root / adaptive_subdir / cid / 'review_candidates.geojson').read_text())['features']
        by_name = {f['properties']['name']: shape(f['geometry']) for f in features}
        indices = [source['source_channel_names'].index(m) for m in MARKERS + ['DAPI', 'CHGA']]
        try:
            pixels = tifffile.memmap(source['source_image'], mode='r')
        except ValueError:
            pixels = tifffile.imread(source['source_image'], key=indices, maxworkers=1)
            indices = list(range(5))
        # Independent context displays use image-wide quantiles; no per-crop scale.
        context_ranges = [np.quantile(pixels[i, ::16, ::16], [.05, .999]) for i in indices[-2:]]
        low = np.asarray(prediction['normalization']['offset'])
        high = low + 12 * np.asarray(prediction['normalization']['scale'])
        fig, axes = plt.subplots(len(names), 4, figsize=(14, 3.5 * len(names)), squeeze=False)
        for row, name in enumerate(names):
            geom = by_name[name]
            a, b, c, d = geom.bounds
            pad = 15 / source['pixel_size_um']
            x0, x1 = max(0, int(a - pad)), min(pixels.shape[2], int(c + pad + 1))
            y0, y1 = max(0, int(b - pad)), min(pixels.shape[1], int(d + pad + 1))
            raw = np.stack([pixels[i, y0:y1, x0:x1] for i in indices[:3]], axis=-1)
            rgb = np.clip((raw - low) / np.maximum(high - low, 1), 0, 1)
            for col in range(4):
                ax = axes[row, col]
                if col < 2:
                    ax.imshow(rgb)
                else:
                    lo, hi = context_ranges[col - 2]
                    ax.imshow(pixels[indices[col + 1], y0:y1, x0:x1], cmap='gray', vmin=lo, vmax=max(hi, lo + 1))
                if col:
                    contour(ax, geom, np.array([x0, y0]))
                ax.set_xlim(0, x1 - x0)
                ax.set_ylim(y1 - y0, 0)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_title([f'{name}: raw hormones', 'Boundary', 'DAPI context', 'CHGA context'][col], fontsize=9)
        fig.suptitle(f'Image reference {cid}: targeted context for uncertain boundaries; no disease category shown', fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .96))
        out = output_root / cid
        out.mkdir(parents=True, exist_ok=True)
        for ext in ['png', 'pdf']:
            fig.savefig(out / f'context_review.{ext}', dpi=130)
        plt.close(fig)
        item = {'case_id': cid, 'candidates': names, 'panel_sha256': sha256(out / 'context_review.png'),
                'context_display_ranges': dict(zip(['DAPI', 'CHGA'], [v.tolist() for v in context_ranges])),
                'status': 'rendered_awaiting_visual_review'}
        write_json(out / 'context_manifest.json', item)
        result.append(item)
    return result


def render_edge_panel(boundaries_root, output_root, cid, edge_width=180):
    """Screen all four field edges for non-biological registration/padding bands."""
    boundaries_root, output_root = Path(boundaries_root), Path(output_root)
    out = output_root / cid
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'edge_review.png').is_file():
        return {'case_id': cid, 'panel': str(out / 'edge_review.png'), 'status': 'already_rendered'}
    source = json.loads((boundaries_root / 'sources' / cid / 'source.json').read_text())
    pred = boundaries_root / 'predictions' / cid / 'prediction_manifest.json'
    if not pred.is_file():
        return {'case_id': cid, 'status': 'prediction_not_ready'}
    prediction = json.loads(pred.read_text())
    indices = [source['source_channel_names'].index(m) for m in MARKERS]
    try:
        image = tifffile.memmap(source['source_image'], mode='r')
    except ValueError:
        image = tifffile.imread(source['source_image'], key=indices, maxworkers=1)
        indices = list(range(3))
    height, width = image.shape[1:]
    lo = np.asarray(prediction['normalization']['offset'])
    hi = lo + 12 * np.asarray(prediction['normalization']['scale'])
    panels = [np.stack([image[i, :edge_width, ::4] for i in indices], -1),
              np.stack([image[i, -edge_width:, ::4] for i in indices], -1),
              np.stack([image[i, ::4, :edge_width].T for i in indices], -1),
              np.stack([image[i, ::4, -edge_width:].T for i in indices], -1)]
    ranges = [[0, width, edge_width, 0], [0, width, height, height-edge_width],
              [0, height, edge_width, 0], [0, height, width, width-edge_width]]
    fig, axes = plt.subplots(4, 1, figsize=(17, 8))
    for ax, data, extent, label in zip(axes, panels, ranges, ['Top', 'Bottom', 'Left (transposed)', 'Right (transposed)']):
        ax.imshow(np.clip((data.astype(float)-lo)/np.maximum(hi-lo, 1), 0, 1), extent=extent, aspect='auto', interpolation='nearest')
        ax.set_title(label, fontsize=10)
        ax.set_ylabel('Source pixel edge axis')
    fig.suptitle(f'Image reference {cid}: four-edge acquisition-artifact screen, INS/GCG/SST\nShort dimension expanded for visibility. Uniform rectangular bands require a technical assessability flag.', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(out / 'edge_review.png', dpi=120)
    plt.close(fig)
    record = {'case_id': cid, 'panel': str(out / 'edge_review.png'), 'panel_sha256': sha256(out / 'edge_review.png'),
              'source_pixel_edge_width': edge_width, 'status': 'rendered_awaiting_visual_review',
              'source_image_modified': False, 'owner_approval_claimed': False}
    write_json(out / 'edge_manifest.json', record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--boundaries-root', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--case-id', action='append', default=[])
    parser.add_argument('--ready-only', action='store_true')
    parser.add_argument('--adaptive-subdir', default='adaptive', help='Explicit geometry generation, e.g. adaptive_v2; output roots must remain separate')
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        raise ValueError('Boundary audit reserves at most four CPU workers')
    if Path(args.adaptive_subdir).name != args.adaptive_subdir:
        raise ValueError('adaptive-subdir must be a direct child directory name')
    donors = args.case_id or sorted(p.name for p in (args.boundaries_root / 'sources').iterdir() if (p / 'source.json').is_file())
    ready = [cid for cid in donors if (args.boundaries_root / args.adaptive_subdir / cid / 'vectorization_manifest.json').is_file()]
    pending = sorted(set(donors) - set(ready))
    if pending and not args.ready_only:
        raise ValueError(f'Boundary outputs are not complete for {pending}')
    args.output_root.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(audit_one, args.boundaries_root, args.output_root, cid, args.adaptive_subdir) for cid in ready]
        records = [future.result() for future in futures]
    write_json(args.output_root / 'audit_manifest.json', {'schema_version': 1,
        'status': 'rendered_awaiting_visual_review' if not pending else 'partial_render_waiting_for_processing',
        'adaptive_subdir': args.adaptive_subdir,
        'donors_requested': donors, 'donors_rendered': ready, 'pending_donors': pending,
        'geometry_passed': [r['case_id'] for r in records if r['geometry_check']['status'] == 'passed'],
        'n_sampled_examples': sum(len(r['selected_examples']) for r in records),
        'visual_review_completed': False, 'owner_approval_claimed': False})


if __name__ == '__main__':
    main()
