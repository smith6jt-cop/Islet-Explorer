#!/usr/bin/env python3
"""Independently check serialized QuPath geometry with Shapely.

QuPath ROI/JTS conversion has finite source-coordinate precision. Report actual
rounding fringes and test their displacement, rather than interpreting tiny
subpixel differences as biologically assessable tissue area.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import box, shape

from prepare_endocrine_candidates import digest, write_json


def overlap_exceeds_rounding(intersection, displacement_px):
    """Two independently rounded edges may meet in a strip up to 2*error wide.

    Reject any intersection that retains interior beyond that strip. The error
    is derived from native precision, never fitted to an observed overlap.
    """
    return not intersection.buffer(-displacement_px).is_empty


def verify(folder):
    manifest = json.loads((folder/'prediction_manifest.json').read_text())
    vector = json.loads((folder/'vectorization_manifest.json').read_text())
    path = folder/'review_candidates.geojson'
    record = next(r for r in vector['outputs'] if r['kind'] == 'boundary_regions')
    if digest(path) != record['sha256']:
        raise ValueError('Changed candidate objects')
    features = json.loads(path.read_text())['features']
    gs = [shape(f['geometry']) for f in features]
    pixel_um = manifest['source_pixel_size_um']
    tree = shapely.STRtree(gs)
    overlap_area, overlap_pairs = 0., 0
    intersections = []
    for i, geom in enumerate(gs):
        for j in tree.query(geom, predicate='intersects'):
            if j >= i:
                continue
            intersection = geom.intersection(gs[j])
            area = intersection.area*pixel_um**2
            if area > 1e-6:
                overlap_pairs += 1
                overlap_area += area
                intersections.append(intersection)
    spec = Path(manifest['source_specification'])
    source_tissue_file = spec.parent/'valid_tissue.geojson'
    if digest(source_tissue_file) != vector['tissue_sha256']:
        raise ValueError('Changed tissue geometry')
    tissue_file = Path(manifest['assessable_tissue_path']) if manifest.get('assessable_tissue_path') else source_tissue_file
    if manifest.get('assessable_tissue_sha256') and digest(tissue_file)!=manifest['assessable_tissue_sha256']:
        raise ValueError('Changed technical assessability geometry')
    tissue = shape(json.loads(tissue_file.read_text()))
    margin = manifest['boundary_config']['edge_margin_um']/pixel_um
    allowed = tissue.intersection(box(margin, margin, manifest['source_width']-margin, manifest['source_height']-margin))
    fringes = [g.difference(allowed) for g in gs]
    outside_area = sum(g.area for g in fringes)*pixel_um**2
    coords = shapely.get_coordinates(fringes)
    max_outside = float(np.max(shapely.distance(shapely.points(coords), allowed), initial=0))*pixel_um
    # Installed QuPath 0.7 GeometryTools uses a 0.01-unit JTS grid and float
    # coordinates. The grid is applied in the physical buffer/clip frame and
    # again in source pixels during ROI conversion. Account for both rounding
    # steps plus one source-coordinate float ULP; verify the native grid in the
    # companion QuPath project/import-gate check.
    if vector['qupath_version'] != '0.7.0':
        raise ValueError('Review native coordinate precision for this QuPath version')
    geometry_grid = .01
    float_ulp_px = float(np.spacing(np.float32(max(manifest['source_width'], manifest['source_height']))))
    displacement_tolerance_um = np.sqrt(2)*.5*geometry_grid*(1+pixel_um) + float_ulp_px*pixel_um
    coordinate_tolerance_px = displacement_tolerance_um/pixel_um
    beyond_precision = sum(overlap_exceeds_rounding(g, coordinate_tolerance_px) for g in intersections)
    areas = np.asarray([g.area*pixel_um**2 for g in gs])
    report = dict(case_id=manifest['case_id'], status='unchecked', objects=len(gs),
        all_geometry_valid=all(g.is_valid for g in gs), positive_overlap_pairs=overlap_pairs,
        positive_overlap_area_um2=overlap_area,
        overlap_pairs_beyond_coordinate_precision=beyond_precision,
        pre_serialization_positive_overlap_pairs=manifest['boundary_summary']['positive_area_overlaps'],
        overlap_policy='Report every >1e-6 um2 intersection; reject interior beyond native coordinate rounding or any pre-serialization overlap',
        native_rounding_area_outside_tissue_and_inner_image_um2=outside_area,
        maximum_outside_displacement_um=max_outside,
        physical_geometry_grid_um=geometry_grid,source_geometry_grid_px=geometry_grid,
        source_coordinate_float_ulp_px=float_ulp_px,
        source_coordinate_tolerance_px=float(coordinate_tolerance_px),
        displacement_tolerance_um=float(displacement_tolerance_um),
        minimum_serialized_area_um2=float(areas.min()) if len(areas) else None,
        serialized_area_um2=float(areas.sum()),
        pre_serialization_area_um2=manifest['boundary_summary']['output_area_um2'],
        review_candidates_sha256=digest(path), original_sources_modified=False)
    if (not report['all_geometry_valid'] or beyond_precision
            or report['pre_serialization_positive_overlap_pairs']
            or max_outside > coordinate_tolerance_px*pixel_um):
        raise ValueError(report)
    if any(g.area*pixel_um**2 + g.length*coordinate_tolerance_px*pixel_um**2 < 50 for g in gs):
        raise ValueError('A fragment is below the minimum beyond coordinate-rounding tolerance')
    if len(gs) != record['objects']:
        raise ValueError('Serialized count differs')
    report['status'] = 'passed'
    write_json(folder/'independent_geometry_check.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folders', nargs='+', type=Path)
    for folder in parser.parse_args().folders:
        print(json.dumps(verify(folder)), flush=True)


if __name__ == '__main__':
    main()
