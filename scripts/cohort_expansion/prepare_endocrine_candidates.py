#!/usr/bin/env python3
"""Prepare a separate, endocrine-seeded boundary review; never validate phenotypes.

Exported hormone thresholds are operational seed thresholds, not reviewed gates.
All seed clusters, including one-cell foci, are retained. Disease labels are unused.
Run with the analysis environment; SAM inference uses the separate SAM environment.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage
import shapely
from shapely.geometry import shape, mapping
from sklearn.cluster import DBSCAN
import tifffile
from convert_images import canonical_channel


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def refuse_withdrawn_candidate(root):
    """A known failed candidate must not be repackaged as a new delivery."""
    path = Path(root)/'quality_status.json'
    if path.exists():
        status = str(json.loads(path.read_text())['status'])
        if status.startswith(('withdrawn', 'failed')):
            raise ValueError(f'Boundary quality withdrawal blocks delivery: {path} ({status})')


def spread_points(points, cap=32):
    """Deterministic spatial coverage; no marker magnitude or disease ordering."""
    if len(points) <= cap:
        return np.arange(len(points))
    chosen = [int(np.argmin(np.sum((points-points.mean(axis=0))**2, axis=1)))]
    distances = np.full(len(points), np.inf)
    while len(chosen) < cap:
        distances = np.minimum(distances, np.sum((points-points[chosen[-1]])**2, axis=1))
        distances[chosen] = -1
        chosen.append(int(np.argmax(distances)))
    return np.array(chosen)


def hormone_points(rendered, centres, radius_px):
    """Move each prompt to nearby hormone signal; centroids may be dark nuclei.

    This proximity rule is an operational prompt, not a replacement cell location.
    Original seed centroids remain the coordinates used for all support counts.
    """
    signal = ndimage.gaussian_filter(rendered.astype(float).max(axis=-1), sigma=1)
    result = []
    for x, y in centres:
        x0, y0 = max(0, int(np.floor(x-radius_px))), max(0, int(np.floor(y-radius_px)))
        x1, y1 = min(signal.shape[1], int(np.ceil(x+radius_px))+1), min(signal.shape[0], int(np.ceil(y+radius_px))+1)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        distance = (xx-x)**2+(yy-y)**2
        local = np.where(distance <= radius_px**2, signal[y0:y1, x0:x1], -np.inf)
        best = np.flatnonzero(local == local.max())
        # For equal signal, choose the nearest pixel; row/column order breaks ties.
        idx = best[np.argmin(distance.ravel()[best])]
        result.append([int(xx.ravel()[idx]), int(yy.ravel()[idx])])
    return np.asarray(result, dtype=float)


def prepare(root, review, out, cid, spacing, eps):
    source = pd.read_csv(root/'source_manifest.csv', dtype=str).set_index('case_id').loc[cid]
    donor = next(d for d in review['donors'] if d['case_id'] == cid)
    exported = json.loads(Path(donor['export_manifest']).read_text())
    if exported['pixel_size_um_x'] != spacing or exported['pixel_size_um_y'] != spacing:
        raise ValueError('Requested spatial calibration differs from the measured source')
    image = Path(donor['source_image']).resolve(strict=True)
    raw_path = root/'raw'/f'{cid}-cells.parquet'
    comp_path = root/'regional'/f'islet_composition_{int(cid):04d}.tsv'
    table = pd.read_csv(comp_path, sep='\t')
    thresholds = {}
    for marker, column in [('INS', 'Ins_thr'), ('GCG', 'Glu_thr'), ('SST', 'Stt_thr')]:
        values = table[column].unique()
        if len(values) != 1 or not np.isfinite(values[0]) or values[0] < 0:
            raise ValueError(f'Expected one finite donor threshold: {cid} {marker}')
        thresholds[marker] = float(values[0])
    cells = pd.read_parquet(raw_path, columns=['object_id', 'X_centroid', 'Y_centroid', 'INS', 'GCG', 'SST'])
    values = cells[list(thresholds)].to_numpy()
    positive = np.isfinite(values).all(axis=1) & (values > np.array(list(thresholds.values()))).any(axis=1)
    with gzip.open(source.geometry, 'rt') as stream:
        source_geometry = json.load(stream)
    tissue, excluded = [], []
    for feature in source_geometry['features']:
        cls = feature['properties'].get('classification', {}).get('name')
        if cls == 'Tissue': tissue.append(shapely.make_valid(shape(feature['geometry'])))
        if cls == 'Exclude': excluded.append(shapely.make_valid(shape(feature['geometry'])))
    if not tissue:
        raise ValueError(f'No tissue mask: {cid}')
    valid = shapely.union_all(tissue).difference(shapely.union_all(excluded))
    xy = cells[['X_centroid', 'Y_centroid']].to_numpy()/spacing
    in_tissue = shapely.intersects_xy(valid, xy[:, 0], xy[:, 1])
    seeds = cells.loc[positive & in_tissue].sort_values('object_id').copy()
    seed_xy = seeds[['X_centroid', 'Y_centroid']].to_numpy()/spacing
    labels = DBSCAN(eps=eps/spacing, min_samples=1).fit_predict(seed_xy)
    groups = sorted((np.flatnonzero(labels == label) for label in np.unique(labels)),
                    key=lambda idx: (float(seed_xy[idx, 1].mean()), float(seed_xy[idx, 0].mean())))
    with tifffile.TiffFile(image) as tf:
        channel_names = [canonical_channel(s) for s in tf.imagej_metadata['Labels']]
    if len(channel_names) != len(set(channel_names)):
        raise ValueError('Ambiguous source channel mapping')
    pixels = tifffile.memmap(image, mode='r')
    if pixels.ndim != 3 or pixels.shape[0] != len(channel_names):
        raise ValueError('Expected a CYX ImageJ TIFF')
    channels = [channel_names.index(m) for m in thresholds]
    display_path = Path(donor['export_manifest']).parent/'sam_prompts.json'
    source_display = json.loads(display_path.read_text())
    display = {d['name'].split(' (')[0]: d for d in source_display['display']}
    levels = [float(display[m]['maxDisplay']) for m in thresholds]
    minima = [float(display[m]['minDisplay']) for m in thresholds]
    if not np.all(np.asarray(levels) > np.asarray(minima)):
        raise ValueError('Invalid fixed display ranges')
    folder = out/cid
    folder.mkdir(parents=True, exist_ok=False)
    (folder/'crops').mkdir()
    prompts = []
    for number, idx in enumerate(groups, 1):
        points = seed_xy[idx]
        name = f'EC_{cid}_{number:04d}'
        seed_record = dict(candidate_id=name, seed_ids=seeds.iloc[idx].object_id.tolist(),
                           seed_xy_px=points.tolist(), seed_count=len(idx))
        if len(idx) >= 3:
            lo, hi = points.min(axis=0)-12/spacing, points.max(axis=0)+12/spacing
            # Local crops retain source resolution unless longer than 900 pixels.
            low = np.maximum(0, np.floor(lo-25/spacing).astype(int))
            high = np.minimum(np.array(pixels.shape[1:][::-1]), np.ceil(hi+25/spacing).astype(int))
            downsample = max(1, int(np.ceil(max(high-low)/900)))
            crop = np.stack([pixels[c, low[1]:high[1]:downsample, low[0]:high[0]:downsample]
                             for c in channels], axis=-1)
            rendered = np.uint8(np.clip((crop-np.asarray(minima))/(np.asarray(levels)-minima), 0, 1)*255)
            Image.fromarray(rendered).save(folder/'crops'/f'{name}.png')
            selection = spread_points(points)
            centres = (points[selection]-low)/downsample
            foreground = hormone_points(rendered, centres, 4/(spacing*downsample))
            seed_record.update(crop=f'crops/{name}.png', offset_px=low.tolist(), downsample=downsample,
                               centroid_points=centres.tolist(), points=foreground.tolist(),
                               box=np.r_[np.maximum(0, (lo-low)/downsample),
                                         np.minimum(np.array(rendered.shape[1::-1])-1, (hi-low)/downsample)].tolist())
        else:
            seed_record['review_class'] = 'Endocrine focus: 1-2 cells'
        prompts.append(seed_record)
    seeds['cluster_id'] = ''
    for prompt in prompts:
        seeds.loc[seeds.object_id.isin(prompt['seed_ids']), 'cluster_id'] = prompt['candidate_id']
    seeds.to_csv(folder/'seed_cells.csv', index=False)
    write_json(folder/'valid_tissue.geojson', mapping(valid))
    spec = dict(case_id=cid, source_image=str(image), image_shape_cyx=list(pixels.shape),
                pixel_size_um=spacing, operational_thresholds=thresholds, threshold_comparator='>',
                thresholds_biologically_validated=False, clustering_eps_um=eps,
                rgb_channels=list(thresholds), source_channel_names=channel_names,
                display_min=minima, display_max=levels,
                display_method='Reuse archived per-image channel ranges; no per-crop contrast fitting',
                source_display=str(display_path), source_display_sha256=digest(display_path),
                prompt_box_padding_um=12, crop_context_um=25, point_cap=32,
                point_placement='Maximum smoothed hormone RGB signal within 4 um of selected seed centroid',
                point_peak_radius_um=4, point_smoothing_sigma_render_pixels=1,
                source_geometry=str(source.geometry), source_geometry_sha256=digest(source.geometry),
                raw_cells=str(raw_path), raw_cells_sha256=digest(raw_path),
                source_thresholds=str(comp_path), source_thresholds_sha256=digest(comp_path),
                raw_cell_count=len(cells), hormone_positive_cells=int(positive.sum()),
                hormone_positive_outside_valid_tissue=int((positive & ~in_tissue).sum()),
                seed_cell_count=len(seeds), cluster_count=len(prompts),
                sam_prompt_count=sum('crop' in p for p in prompts), prompts=prompts,
                acceptance='Unreviewed candidates only; not phenotype validation or analytical acceptance')
    write_json(folder/'prompts.json', spec)
    print(f'{cid}: {len(seeds)} seed cells; {spec["sam_prompt_count"]} SAM prompts; '
          f'{len(prompts)-spec["sam_prompt_count"]} small foci retained', flush=True)


def reprompt(source, destination, root, review, cid, spacing, eps):
    """Reuse independently copied, unchanged renders when only prompt placement changes."""
    spec = json.loads((source/'prompts.json').read_text())
    donor = next(d for d in review['donors'] if d['case_id'] == cid)
    if (spec['case_id'] != cid or spec['pixel_size_um'] != spacing or spec['clustering_eps_um'] != eps
        or Path(spec['raw_cells']).resolve() != root/'raw'/f'{cid}-cells.parquet'
        or Path(spec['source_image']).resolve() != Path(donor['source_image']).resolve()):
        raise ValueError('Cannot reuse renders from different inputs or clustering')
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copytree(source/'crops', destination/'crops')
    for filename in ['seed_cells.csv', 'valid_tissue.geojson']:
        shutil.copy2(source/filename, destination/filename)
    for prompt in spec['prompts']:
        if 'crop' not in prompt: continue
        rendered = np.asarray(Image.open(destination/prompt['crop']).convert('RGB'))
        centres = np.asarray(prompt.get('centroid_points', prompt['points']))
        prompt['centroid_points'] = centres.tolist()
        prompt['points'] = hormone_points(rendered, centres, 4/(spacing*prompt['downsample'])).tolist()
    spec.update(point_placement='Maximum smoothed hormone RGB signal within 4 um of selected seed centroid',
                point_peak_radius_um=4, point_smoothing_sigma_render_pixels=1,
                reused_render_source=str(source), source_prompts_sha256=digest(source/'prompts.json'))
    write_json(destination/'prompts.json', spec)
    print(f'{cid}: reprompted {spec["sam_prompt_count"]} independently copied crops', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root', required=True, type=Path)
    p.add_argument('--review-manifest', required=True, type=Path)
    p.add_argument('--output-root', required=True, type=Path)
    p.add_argument('--case-ids', nargs='+', default=['78'])
    p.add_argument('--pixel-size-um', type=float, default=.3774)
    p.add_argument('--cluster-distance-um', type=float, default=20)
    p.add_argument('--reuse-candidate-root', type=Path)
    a = p.parse_args()
    if a.pixel_size_um <= 0 or a.cluster_distance_um <= 0:
        raise ValueError('Spatial units must be positive')
    review = json.loads(a.review_manifest.read_text())
    for cid in a.case_ids:
        if a.reuse_candidate_root:
            reprompt(a.reuse_candidate_root.resolve(strict=True)/cid, a.output_root.resolve()/cid,
                     a.input_root.resolve(strict=True), review, cid, a.pixel_size_um, a.cluster_distance_um)
        else:
            prepare(a.input_root.resolve(strict=True), review, a.output_root.resolve(), cid,
                    a.pixel_size_um, a.cluster_distance_um)


if __name__ == '__main__': main()
