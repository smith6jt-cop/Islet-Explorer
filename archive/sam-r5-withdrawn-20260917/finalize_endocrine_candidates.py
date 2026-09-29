#!/usr/bin/env python3
"""Make QuPath review overlays, accounting for every hormone seed and old outline."""
import argparse
import gzip
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape, mapping, MultiPoint

# Retain explicit historical execution from this archive with the shared
# preparation helpers that remain in scripts/cohort_expansion.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts' / 'cohort_expansion'))
from prepare_endocrine_candidates import digest, write_json, refuse_withdrawn_candidate


def finish(folder):
    refuse_withdrawn_candidate(folder.parent)
    if (folder/'review_candidates.geojson').exists():
        raise ValueError('Do not overwrite review candidates')
    spec = json.loads((folder/'prompts.json').read_text())
    cid, scale = spec['case_id'], spec['pixel_size_um']
    audit = json.loads((folder/'prediction_audit.json').read_text())
    if (audit['case_id'] != cid or folder.name != cid
        or audit['prompts_sha256'] != digest(folder/'prompts.json')
        or audit['output_sha256'] != digest(folder/'sam_predictions.geojson')):
        raise ValueError('Prediction identity or prompt provenance differs')
    if digest(spec['raw_cells']) != spec['raw_cells_sha256']:
        raise ValueError('Source cells changed since prompt construction')
    if digest(spec['source_geometry']) != spec['source_geometry_sha256']:
        raise ValueError('Source geometry changed since prompt construction')
    predictions = json.loads((folder/'sam_predictions.geojson').read_text())['features']
    geometries = np.array([shape(p['geometry']) for p in predictions], dtype=object)
    tree = shapely.STRtree(geometries)
    seeds = pd.read_csv(folder/'seed_cells.csv')
    seed_xy = seeds[['X_centroid', 'Y_centroid']].to_numpy()/scale
    points = shapely.points(seed_xy)
    seed_tree = shapely.STRtree(points)
    cells = pd.read_parquet(spec['raw_cells'], columns=['object_id', 'X_centroid', 'Y_centroid', 'INS', 'GCG', 'SST'])
    cell_points = shapely.points(cells[['X_centroid', 'Y_centroid']].to_numpy()/scale)
    cell_tree = shapely.STRtree(cell_points)
    # Union intersecting predictions; do not introduce dilation, convex hulls, or bridges.
    parts = [g for g in shapely.get_parts(shapely.union_all(geometries)) if g.geom_type == 'Polygon']
    parts.sort(key=lambda g: (g.centroid.y, g.centroid.x))
    features, rows, represented = [], [], np.zeros(len(seeds), dtype=bool)
    for i, geom in enumerate(parts, 1):
        if not geom.is_valid or geom.area <= 0:
            raise ValueError('Invalid final geometry')
        selected = seed_tree.query(geom, predicate='intersects')
        source_ids = tree.query(geom, predicate='intersects')
        flags = [predictions[j]['properties'] for j in source_ids]
        edge = any(p['crop_edge_contact'] for p in flags)
        recall = min(p['prompted_seed_recall'] for p in flags)
        inverse = any(p.get('inverse_or_flat_signal_contrast', False) for p in flags)
        if len(selected) == 0:
            # Retain in the prediction file/audit; never create a seedless review islet.
            continue
        represented[selected] = True
        cell_ids = cell_tree.query(geom, predicate='intersects')
        cls = ('Boundary review: sparse support' if len(selected) < 3 else
               'Boundary review: edge, partial or low contrast' if edge or recall < .8 or inverse else
               'Islet candidate: endocrine supported')
        name = f'New_{cid}_{i:04d}'
        measurements = {'Candidate: hormone-positive seed cells': len(selected),
                        'Candidate: total cells in new boundary': len(cell_ids),
                        'Candidate: area um2': geom.area*scale**2,
                        'Candidate: source SAM prompts': len(source_ids),
                        'Candidate: minimum prompted seed recall': recall,
                        'Candidate: crop edge contact': int(edge),
                        'Candidate: inverse or flat hormone contrast': int(inverse),
                        'Islet: identity reviewed': 0}
        features.append(dict(type='Feature', geometry=mapping(geom), properties=dict(
            objectType='annotation', name=name,
            classification=dict(name=cls, color=[0, 210, 120] if cls.startswith('Islet candidate') else [240, 190, 0]),
            measurements=measurements, source_prompt_ids=[p['candidate_id'] for p in flags],
            review_status='unreviewed; operational hormone thresholds are not validated phenotype gates')))
        rows.append(dict(candidate_id=name, classification=cls, area_um2=geom.area*scale**2,
                         hormone_seed_cells=len(selected), total_cells=len(cell_ids),
                         source_prompt_count=len(source_ids), crop_edge_contact=edge,
                         inverse_or_flat_hormone_contrast=inverse,
                         prompted_seed_recall=recall, decision='unreviewed', review_note=''))
    focus_count = 0
    for prompt in spec['prompts']:
        idx = np.flatnonzero(seeds.object_id.isin(prompt['seed_ids']).to_numpy() & ~represented)
        if not len(idx): continue
        is_small = prompt['seed_count'] < 3
        cls = 'Endocrine focus: 1-2 cells' if is_small else 'Endocrine seed: boundary unresolved'
        features.append(dict(type='Feature', geometry=mapping(MultiPoint(seed_xy[idx])), properties=dict(
            objectType='annotation', name='Focus_'+prompt['candidate_id'],
            classification=dict(name=cls, color=[70, 160, 255]),
            measurements={'Candidate: hormone-positive seed cells': len(idx), 'Islet: identity reviewed': 0},
            review_status='Point markers only, not islet boundaries')))
        focus_count += 1
        represented[idx] = True
    if not represented.all(): raise ValueError('Some endocrine seeds disappeared from the review')
    write_json(folder/'review_candidates.geojson', dict(type='FeatureCollection', features=features))
    pd.DataFrame(rows).to_csv(folder/'candidate_review.csv', index=False)
    with gzip.open(spec['source_geometry'], 'rt') as stream:
        originals = json.load(stream)['features']
    old = [f for f in originals if f['properties'].get('classification', {}).get('name') == 'Islet']
    boundary_features = [f for f in features if f['geometry']['type'] in ['Polygon', 'MultiPolygon']]
    union = shapely.union_all([shape(f['geometry']) for f in boundary_features])
    comparison = []
    for original in old:
        geom = shape(original['geometry'])
        comparison.append(dict(previous_islet=original['properties']['name'],
            old_hormone_seed_cells=len(seed_tree.query(geom, predicate='intersects')),
            old_area_um2=geom.area*scale**2,
            fraction_covered_by_new_boundaries=geom.intersection(union).area/geom.area if geom.area else None))
    pd.DataFrame(comparison).to_csv(folder/'previous_islet_comparison.csv', index=False)
    summary = dict(case_id=cid, status='review candidates; biological acceptance pending',
        old_islet_count=len(old), sam_prompt_count=spec['sam_prompt_count'],
        new_boundary_count=len(rows), endocrine_supported_candidates=sum(r['classification'].startswith('Islet candidate') for r in rows),
        boundaries_needing_extra_review=sum(not r['classification'].startswith('Islet candidate') for r in rows),
        separate_point_focus_annotations=focus_count, seed_cell_count=len(seeds),
        every_seed_represented=bool(represented.all()),
        no_seedless_boundary=all(r['hormone_seed_cells'] > 0 for r in rows),
        geometry_units='source image pixels', pixel_size_um=scale,
        overlapping_masks='union only; disconnected components retained; no morphological bridging',
        small_foci='Point annotations retain unmasked seed clusters, including single cells',
        old_objects_preserved=True, accepted_islets=0,
        analysis_memberships_rebuilt=False, analytical_acceptance=False,
        overlay_sha256=digest(folder/'review_candidates.geojson'))
    write_json(folder/'candidate_summary.json', summary)
    print(json.dumps(summary), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--case-ids', nargs='+', default=['78'])
    a = p.parse_args()
    for cid in a.case_ids: finish(a.root/cid)


if __name__ == '__main__': main()
