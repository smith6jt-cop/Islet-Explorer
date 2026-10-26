#!/usr/bin/env python3
"""Analysis-only edge eligibility from original cell polygons and saved pre-clip masks.

Never modifies the accepted app or source run. The retained classifier mask is
read BEFORE script 1's 5 um inset and morphology. Whole connected components are
excluded conservatively: every final object intersecting a border component is
withheld, including interior descendants after a split. Eight-connectivity and
one classifier-pixel matching tolerance deliberately err toward exclusion.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import ndimage
import shapely
from shapely.geometry import box, shape
import tifffile

from h1_features_20261022 import sha256, write_json

VERSION = 'h1-edge-20261024-v3'
CELL_TOLERANCE_PX = 1.0
INSET_UM = 5.0
RADII = (20, 100, 250)


def border_contact(geom, width, height, tolerance=CELL_TOLERANCE_PX):
    """Complete original polygon bounds, not its centroid or a clipped raster."""
    if geom is None or geom.is_empty:
        return None
    x0, y0, x1, y1 = geom.bounds
    if not np.isfinite([x0, y0, x1, y1]).all():
        return None
    return bool(x0 <= tolerance or y0 <= tolerance
                or x1 >= width - tolerance or y1 >= height - tolerance)


def restrict_to_tissue(mask, tissue, ds):
    """Match the original classifier's selected Tissue before component labeling.

    Keep one classifier pixel of boundary tolerance for raster support. Tissue
    boundaries themselves are not acquisition edges.
    """
    result = mask.copy()
    support = tissue.buffer(ds)
    shapely.prepare(support)
    for y in range(0, mask.shape[0], 128):
        iy, ix = np.nonzero(result[y:y + 128])
        keep = shapely.intersects_xy(support, (ix + .5) * ds, (iy + y + .5) * ds)
        result[iy[~keep] + y, ix[~keep]] = False
    return result


def damaged_member_mask(cells):
    """Propagate observed border contact, including source-QC-rejected cells.

    Missing geometry in an already QC-rejected segmentation is not evidence
    of image truncation. Missing geometry in an otherwise eligible cell leaves
    its object unassessable and therefore withheld.
    """
    contact = cells[['cell_edge', 'nucleus_edge']].eq(True).any(axis=1)
    return contact | (cells.geometry_missing & cells.qc_analysis_eligible)


def label_preclip(mask, ds, width, height, px, tolerance=CELL_TOLERANCE_PX):
    labels, n = ndimage.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    # Include any classifier pixel whose square reaches the inset plus the
    # physical-border tolerance. Do not infer damage from the final outline.
    pad = INSET_UM / px + tolerance
    left = np.arange(mask.shape[1]) * ds < pad
    right = (np.arange(mask.shape[1]) + 1) * ds > width - pad
    top = np.arange(mask.shape[0]) * ds < pad
    bottom = (np.arange(mask.shape[0]) + 1) * ds > height - pad
    ids = np.unique(np.concatenate([labels[:, left | right].ravel(), labels[top | bottom].ravel()]))
    flags = np.zeros(n + 1, dtype=bool)
    flags[ids] = True
    flags[0] = False
    return labels, flags


def component_ids(geom, labels, ds):
    """Conservative raster-component intersection, including pixel support.

    Buffer by one grid pixel to account for vector/raster boundary differences.
    This can associate extra ancestors, which only withholds extra objects.
    """
    g = geom.buffer(ds)
    x0, y0, x1, y1 = g.bounds
    h, w = labels.shape
    xa, xb = max(0, int(np.floor(x0 / ds))), min(w, int(np.ceil(x1 / ds)))
    ya, yb = max(0, int(np.floor(y0 / ds))), min(h, int(np.ceil(y1 / ds)))
    if xa >= xb or ya >= yb:
        return np.array([], dtype=int)
    local = labels[ya:yb, xa:xb]
    iy, ix = np.nonzero(local)
    hit = shapely.intersects_xy(g, (ix + xa + .5) * ds, (iy + ya + .5) * ds)
    return np.unique(local[iy[hit], ix[hit]])


def neighborhood_support(core_um, acquisition_um, radius):
    requested = core_um.buffer(radius)
    observed = requested.intersection(acquisition_um)
    missing = requested.difference(acquisition_um).area
    # Floating-point area tolerance only; no partial image support is accepted.
    complete = missing <= max(1e-7, requested.area * 1e-12)
    return bool(complete), float(requested.area), float(observed.area)


def image_metadata(run, entry):
    path = run / 'qupath' / 'data' / str(entry['entry_id']) / 'server.json'
    md = json.loads(path.read_text())['metadata']
    width, height = int(md['width']), int(md['height'])
    cal = md['pixelCalibration']
    px, py = float(cal['pixelWidth']['value']), float(cal['pixelHeight']['value'])
    if not np.isclose(px, .3774) or not np.isclose(py, px):
        raise ValueError('Unexpected calibration: refuse the legacy isotropic feature extractor')
    ome_calibration = {}
    with tifffile.TiffFile(entry['image_path']) as tif:
        series = tif.series[0]
        dims = dict(zip(series.axes, series.shape))
        if (dims['X'], dims['Y']) != (width, height):
            raise ValueError(f"{entry['image_name']}: TIFF/QuPath dimensions differ")
        ome = tif.ome_metadata
        if ome:
            import xml.etree.ElementTree as ET
            pixels = next(ET.fromstring(ome).iter('{http://www.openmicroscopy.org/Schemas/OME/2016-06}Pixels'), None)
            if pixels is not None:
                for axis, expected in [('X', px), ('Y', py)]:
                    value = pixels.attrib.get('PhysicalSize' + axis)
                    unit = pixels.attrib.get('PhysicalSize' + axis + 'Unit', 'µm')
                    if value:
                        factors = {'µm': 1., 'um': 1., 'mm': 1000., 'nm': .001}
                        physical_um = float(value) * factors.get(unit, np.nan)
                        ome_calibration[axis] = physical_um
                        # The accepted pipeline rounds physical calibration to
                        # four decimal places; retain that coordinate contract.
                        if not np.isfinite(physical_um) or round(physical_um, 4) != round(expected, 4):
                            raise ValueError(f'OME calibration mismatch: {axis} {value} {unit}')
    return width, height, px, path, dict(tiff_axes=series.axes, tiff_shape=list(series.shape),
                                        ome_metadata_present=bool(ome), ome_calibration_um=ome_calibration, pixel_um=px)


def process(args):
    run_string, app_string, out_string, entry = args
    run, app, out = map(Path, (run_string, app_string, out_string))
    tag, case = entry['image_name'], int(entry['case_id'])
    if case == 78:
        raise ValueError('Excluded donor 78')
    start = time.monotonic()
    paths = [run / 'cell_geojson' / f'{tag}.geojson', run / 'islets' / f'{tag}.geojson',
             run / 'typing/context' / f'{tag}.parquet', run / 'fit' / f'{tag}.mask.npy',
             run / 'fit' / f'{tag}.json', run / 'planes' / f'{tag}.planes.json',
             app / 'cell_types' / f'{case}.parquet', run / 'provenance' / tag / 'fit.json',
             run / 'qupath/classifiers/pixel_classifiers' / f'Islets_{tag}.json']
    width, height, px, server_path, meta = image_metadata(run, entry)
    paths.append(server_path)
    tissue_path = run / 'inputs/tissue' / f'{tag}.geojson'
    paths.append(tissue_path)
    fit = json.loads(paths[4].read_text())
    receipt = json.loads(paths[7].read_text())
    if receipt['status'] != 'pass':
        raise ValueError(f'{tag}: classifier fit did not pass')
    for p in (paths[4], paths[8]):
        if sha256(p) != receipt['outputs'][str(p.resolve())]:
            raise ValueError(f'{tag}: classifier/report no longer matches source receipt')
    plane = json.loads(paths[5].read_text())
    ds = float(plane['downsample'])
    mask = np.load(paths[3])
    if mask.dtype != bool or mask.shape != (plane['height'], plane['width']):
        raise ValueError(f'{tag}: invalid retained pre-clip mask')
    if not np.isclose(ds * px, fit['resolution_um']) or abs(mask.shape[1] * ds - width) > ds or abs(mask.shape[0] * ds - height) > ds:
        raise ValueError(f'{tag}: classifier mask transform does not match image')
    tissue_features = json.loads(tissue_path.read_text())['features']
    tissue = shapely.union_all([shapely.make_valid(shape(f['geometry'])) for f in tissue_features
                               if (f.get('properties', {}).get('classification') or {}).get('name') == 'Tissue'])
    if tissue.is_empty:
        raise ValueError(f'{tag}: missing original Tissue selection')
    original_foreground = int(mask.sum())
    mask = restrict_to_tissue(mask, tissue, ds)
    tissue_foreground = int(mask.sum())
    labels, component_edge = label_preclip(mask, ds, width, height, px)
    del mask
    source = json.loads(paths[0].read_text())
    records = []
    for f in source['features']:
        g = shape(f['geometry']) if f.get('geometry') else None
        ng = shape(f['nucleusGeometry']) if f.get('nucleusGeometry') else None
        cell, nucleus = border_contact(g, width, height), border_contact(ng, width, height)
        missing = cell is None or nucleus is None
        reasons = [r for r, condition in [('cell_image_border', cell), ('nucleus_image_border', nucleus),
                                         ('missing_geometry', missing)] if condition]
        records.append(dict(object_id=f['id'], case_id=case, cell_edge=cell, nucleus_edge=nucleus,
                            geometry_missing=missing, edge_analysis_eligible=not reasons,
                            exclusion_reason=';'.join(reasons)))
    del source
    cells = pd.DataFrame(records)
    if cells.object_id.duplicated().any():
        raise ValueError(f'{tag}: duplicate source UUID')
    context = pd.read_parquet(paths[2], columns=['object_id', 'islet_num', 'endocrine_object_id', 'qc_analysis_eligible'])
    cells = cells.merge(context, on='object_id', how='left', validate='one_to_one')
    if cells.qc_analysis_eligible.isna().any():
        raise ValueError(f'{tag}: incomplete source QC join')
    # Source context includes QC-excluded cells too; do not test object damage
    # using only the subset which survived typing and assembly.
    served = pd.read_parquet(paths[6], columns=['object_id', 'phenotype'])
    if served.object_id.duplicated().any() or not set(served.object_id).issubset(set(cells.object_id)):
        raise ValueError(f'{tag}: served cells do not join uniquely to original segmentation')
    cells['served'] = cells.object_id.isin(served.object_id)
    cells = cells.merge(served, on='object_id', how='left', validate='one_to_one')
    bad_members = set(cells.loc[damaged_member_mask(cells), 'endocrine_object_id'].dropna().astype(str))
    objs = json.loads(paths[1].read_text())['features']
    acquisition = box(0, 0, width * px, height * px)
    rows = []
    for f in objs:
        p = f.get('properties') or {}
        cls = (p.get('classification') or {}).get('name')
        if cls not in ('Islet', 'EndocrineCluster', 'EndocrineSingle'):
            continue
        g = shape(f['geometry'])
        if not g.is_valid:
            g = shapely.make_valid(g)
        ancestors = component_ids(g, labels, ds)
        preclip = bool(component_edge[ancestors].any())
        unknown = not len(ancestors)
        member = str(f['id']) in bad_members
        direct = border_contact(g, width, height)
        reasons = [r for r, flag in [('preclip_border_component', preclip), ('edge_member_cell', member),
                                    ('final_outline_image_border', direct), ('preclip_ancestry_unresolved', unknown)] if flag]
        complete = not reasons
        row = dict(case_id=case, image=tag, object_id=f['id'], object_name=p['name'], object_class=cls,
                   core_complete=complete, exclusion_reason=';'.join(reasons),
                   preclip_edge_component=preclip, edge_member_cell=member,
                   preclip_components=';'.join(map(str, ancestors)), preclip_ancestry_unresolved=unknown,
                   final_outline_image_border=direct)
        from shapely.affinity import scale
        um = scale(g, xfact=px, yfact=px, origin=(0, 0))
        row['core_area_um2'] = um.area
        row['distance_to_image_edge_um'] = um.distance(acquisition.boundary)
        for radius in RADII:
            support, requested, observed = neighborhood_support(um, acquisition, radius)
            row[f'neighborhood_complete_{radius}'] = support
            row[f'neighborhood_requested_area_{radius}_um2'] = requested
            row[f'neighborhood_observed_area_{radius}_um2'] = observed
        rows.append(row)
    objects = pd.DataFrame(rows)
    if objects.object_name.duplicated().any():
        raise ValueError(f'{tag}: duplicate object name')
    cells.to_parquet(out / 'cells' / f'{case}.parquet', index=False)
    objects.to_parquet(out / 'objects' / f'{case}.parquet', index=False)
    audit = dict(version=VERSION, case_id=case, image=tag, image_width=width, image_height=height, **meta,
                 source_cells=len(cells), served_cells=len(served),
                 source_edge_cells=int((~cells.edge_analysis_eligible).sum()),
                 source_border_contact_cells=int(cells[['cell_edge', 'nucleus_edge']].eq(True).any(axis=1).sum()),
                 source_qc_rejected_missing_geometry=int((cells.geometry_missing & ~cells.qc_analysis_eligible).sum()),
                 served_edge_cells=int((cells.served & ~cells.edge_analysis_eligible).sum()),
                 source_objects=len(objects), excluded_objects=int((~objects.core_complete).sum()),
                 unresolved_ancestry=int(objects.preclip_ancestry_unresolved.sum()),
                 preclip_edge_components=int(component_edge.sum()),
                 preclip_foreground_pixels=original_foreground, tissue_selected_foreground_pixels=tissue_foreground,
                 source_sha256={str(p.resolve()): sha256(p) for p in paths},
                 image_source=dict(path=entry['image_path'], bytes=Path(entry['image_path']).stat().st_size,
                                   mtime_ns=Path(entry['image_path']).stat().st_mtime_ns),
                 seconds=round(time.monotonic() - start, 2))
    write_json(out / 'audits' / f'{case}.json', audit)
    print(f"{tag}: {audit['served_edge_cells']} served edge cells; {audit['excluded_objects']}/{len(objects)} objects withheld; {audit['unresolved_ancestry']} unknown ancestors", flush=True)
    return audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--app-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--cases', type=int, nargs='+')
    a = p.parse_args()
    run, app, out = a.source_root.resolve(), a.app_root.resolve(), a.output.resolve()
    for root in (run, app):
        if out == root or out.is_relative_to(root) or root.is_relative_to(out):
            raise ValueError('Output must be separate from source/app data')
    for sub in ('cells', 'objects', 'audits'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    entries = json.loads((run / 'project_manifest.json').read_text())['entries']
    if len(entries) != 20 or any(int(e['case_id']) == 78 for e in entries):
        raise ValueError('Not the accepted 20-donor cohort')
    if a.cases:
        entries = [e for e in entries if int(e['case_id']) in a.cases]
    with ProcessPoolExecutor(a.workers) as pool:
        audits = list(pool.map(process, [(str(run), str(app), str(out), e) for e in entries]))
    objects = pd.concat([pd.read_parquet(out / 'objects' / f"{r['case_id']}.parquet") for r in audits], ignore_index=True)
    objects.to_csv(out / 'object_eligibility.csv', index=False)
    pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, (dict, list))} for r in audits]).to_csv(out / 'donor_exclusions.csv', index=False)
    write_json(out / 'edge_manifest.json', dict(version=VERSION, argv=sys.argv, script_sha256=sha256(__file__),
        donors=len(audits), cell_tolerance_source_pixels=CELL_TOLERANCE_PX, inset_um=INSET_UM,
        object_rule='restrict retained classifier mask to original Tissue plus one classifier-pixel support tolerance; exclude complete pre-inset components reaching inset+one source pixel; propagate by conservative one-classifier-pixel overlap; exclude border-member cells and unresolved ancestors',
        component_connectivity=8, classifier_matching_tolerance_pixels=1,
        acquisition_support='TIFF/QuPath full image rectangle; no independently recorded internal acquisition-gap mask available; internal gaps are not inferred from staining',
        input_manifest_sha256=sha256(run / 'project_manifest.json'), source_run=str(run), app_root=str(app),
        audits=[str(out / 'audits' / f"{r['case_id']}.json") for r in audits]))


if __name__ == '__main__':
    main()
