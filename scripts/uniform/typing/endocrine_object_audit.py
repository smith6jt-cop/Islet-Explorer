#!/usr/bin/env python3
"""Audit endocrine objects without conflating classifier regions and cell means.

Read the selected levels, REDSEA measurements, context membership and final
object inventory directly, before identity typing. Positivity is the same
strict log1p comparison used by identity.py. At least one analysis-eligible
member must be positive for INS, GCG, SST or CHGA; other member cells may be
negative. Empty objects and EndocrineSingle objects without exactly one
analysis-eligible member are separate membership failures.

This is a consistency gate, not a replay of the pixel classifier. Failure does
not prove absence of raw image signal and does not change levels or labels.
Reports must be written to a fresh directory.

An explicit classifier-lineage mode instead validates the preserved positive
pixel-object lineage and Single whole-cell association. It requires both pinned
historical evidence receipts and checks their actual input bytes. The original
corrected-cell test remains unchanged in the report as a diagnostic. Neither
mode changes a cell identity, centroid membership, geometry or measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MARKERS = ('INS', 'GCG', 'SST', 'CHGA')
CLASSES = ('Islet', 'EndocrineCluster', 'EndocrineSingle')
SCHEMA = 'uniform.endocrine_object_audit/2'


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def object_inventory(path: Path) -> pd.DataFrame:
    raw = json.loads(path.read_text())
    features = raw['features'] if isinstance(raw, dict) else raw
    rows = []
    for feature in features:
        props = feature.get('properties') or {}
        label = (props.get('classification') or {}).get('name')
        if label not in CLASSES:
            continue
        if not feature.get('geometry'):
            raise ValueError(f'{path}: endocrine object has no geometry')
        rows.append({'endocrine_object_id': feature.get('id'),
                     'endocrine_object_class': label, 'object_name': props.get('name')})
    frame = pd.DataFrame(rows, columns=['endocrine_object_id', 'endocrine_object_class', 'object_name'])
    if (frame.endocrine_object_id.isna().any() or frame.endocrine_object_id.eq('').any()
            or frame.endocrine_object_id.duplicated().any()):
        raise ValueError(f'{path}: missing or duplicate endocrine object IDs')
    return frame


def features(path: Path) -> list[dict]:
    raw = json.loads(path.read_text())
    return raw['features'] if isinstance(raw, dict) else raw


def streamed_cell_features(path: Path):
    """Read the frozen QuPath FeatureCollection one feature at a time.

    The exported header is explicit; unknown layouts fail rather than silently
    dropping records. No geometry is repaired or rounded by this reader.
    """
    decoder = json.JSONDecoder()
    with path.open() as stream:
        buffer = stream.read(1 << 20)
        header = re.match(r'\s*\{\s*"type"\s*:\s*"FeatureCollection"\s*,\s*"features"\s*:\s*\[', buffer)
        if header is None:
            raise ValueError(f'{path}: unsupported cell FeatureCollection header')
        position = header.end()
        after_feature = False
        while True:
            while True:
                while position < len(buffer) and buffer[position].isspace():
                    position += 1
                if position < len(buffer):
                    break
                buffer, position = stream.read(1 << 20), 0
                if not buffer:
                    raise ValueError(f'{path}: truncated cell FeatureCollection')
            if buffer[position] == ']':
                tail = buffer[position + 1:] + stream.read()
                if tail.strip() != '}':
                    raise ValueError(f'{path}: unexpected FeatureCollection trailer')
                return
            if after_feature:
                if buffer[position] != ',':
                    raise ValueError(f'{path}: missing feature separator')
                position += 1
            while True:
                while position < len(buffer) and buffer[position].isspace():
                    position += 1
                try:
                    feature, position = decoder.raw_decode(buffer, position)
                    break
                except json.JSONDecodeError:
                    more = stream.read(1 << 20)
                    if not more:
                        raise ValueError(f'{path}: malformed or truncated cell feature')
                    buffer = buffer[position:] + more
                    position = 0
            if not isinstance(feature, dict) or feature.get('type') != 'Feature':
                raise ValueError(f'{path}: invalid cell feature')
            yield feature
            after_feature = True
            if position >= 1 << 20:
                buffer, position = buffer[position:], 0


def geometry_bounds(geometry: dict) -> tuple:
    if geometry.get('type') not in ('Polygon', 'MultiPolygon'):
        raise ValueError('Whole-cell association requires polygonal geometry')
    polygons = geometry['coordinates'] if geometry['type'] == 'MultiPolygon' else [geometry['coordinates']]
    points = np.asarray([point[:2] for polygon in polygons for ring in polygon for point in ring], dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not len(points) or not np.isfinite(points).all():
        raise ValueError('Empty or nonfinite polygon coordinates')
    return tuple(np.concatenate([points.min(axis=0), points.max(axis=0)]))


def original_mask_support(geometry: dict, mask: np.ndarray) -> int:
    """Classifier-pixel centers, on the preserved 2x-downsampled image grid."""
    from skimage.draw import polygon as raster_polygon
    polygons = geometry['coordinates'] if geometry['type'] == 'MultiPolygon' else [geometry['coordinates']]
    bounds = np.asarray(geometry_bounds(geometry)) / 2 - .5
    x0, y0 = np.maximum(np.floor(bounds[:2]).astype(int), 0)
    x1, y1 = np.minimum(np.ceil(bounds[2:]).astype(int) + 1, mask.shape[::-1])
    if x1 <= x0 or y1 <= y0:
        return 0
    inside = np.zeros((y1 - y0, x1 - x0), dtype=bool)
    for component in polygons:
        part = np.zeros_like(inside)
        for i, ring in enumerate(component):
            xy = np.asarray(ring, dtype=float)[:, :2] / 2 - .5 - (x0, y0)
            rr, cc = raster_polygon(xy[:, 1], xy[:, 0], shape=part.shape)
            part[rr, cc] = i == 0
        inside |= part
    return int(np.count_nonzero(mask[y0:y1, x0:x1][inside]))


def load_lineage_evidence(classifier_run: Path, classifier_replay: Path, classifier_replay_sha256: str,
                          geometry_reference: Path, geometry_provenance: Path,
                          geometry_provenance_sha256: str, tags: list[str]) -> dict:
    """Validate the two explicit historical pins before considering their claims."""
    for path, expected in ((classifier_replay, classifier_replay_sha256),
                           (geometry_provenance, geometry_provenance_sha256)):
        if not re.fullmatch(r'[0-9a-f]{64}', expected) or sha256(path) != expected:
            raise ValueError(f'{path}: lineage evidence does not match its explicit SHA256 pin')
    replay = json.loads(classifier_replay.read_text())
    provenance = json.loads(geometry_provenance.read_text())
    replay_rows = replay.get('images', [])
    geometry_rows = provenance.get('objects', [])
    def keyed(rows, description):
        if not isinstance(rows, list) or not rows or any(not isinstance(row, dict) for row in rows):
            raise ValueError(f'{description}: invalid image inventory')
        result = {row['image']: row for row in rows}
        if len(result) != len(rows):
            raise ValueError(f'{description}: duplicate image inventory')
        return result
    rp, gp = keyed(replay_rows, 'classifier replay'), keyed(geometry_rows, 'geometry provenance')
    if len(rp) != 20 or set(rp) != set(gp) or not set(tags) <= set(rp):
        raise ValueError('Lineage receipts must cover the same complete 20-image cohort and every audited image')
    if any(row.get('classifier_replay') != 'exact' for row in rp.values()):
        raise ValueError('Historical classifier replay is not exact for every image')
    return {'classifier_run': classifier_run.resolve(), 'geometry_reference': geometry_reference.resolve(),
            'replay': rp, 'geometry': gp,
            'classifier_replay': {'path': str(classifier_replay.resolve()), 'sha256': classifier_replay_sha256},
            'geometry_provenance': {'path': str(geometry_provenance.resolve()), 'sha256': geometry_provenance_sha256}}


def verify_recorded_source(path: Path, suffix: str, records: dict, sources: dict) -> str:
    matches = [digest for name, digest in records.items() if str(name).endswith('/' + suffix)]
    if len(matches) != 1 or not re.fullmatch(r'[0-9a-f]{64}', matches[0]):
        raise ValueError(f'{path}: missing or ambiguous original classifier source hash')
    actual = sha256(path)
    if actual != matches[0]:
        raise ValueError(f'{path}: original classifier input changed')
    sources[str(path)] = actual
    return actual


def audit_classifier_lineage(run: Path, tag: str, evidence: dict, rs: pd.DataFrame,
                              ctx: pd.DataFrame, sources: dict) -> tuple[pd.DataFrame, dict]:
    import shapely
    from shapely.geometry import shape
    base, reference = evidence['classifier_run'], evidence['geometry_reference']
    recorded = evidence['replay'][tag].get('source_sha256', {})
    paths = {'fit': base / 'fit' / f'{tag}.json', 'mask': base / 'fit' / f'{tag}.mask.npy',
             'planes': base / 'planes' / f'{tag}.planes.tif',
             'plane_metadata': base / 'planes' / f'{tag}.planes.json',
             'original_outlines': base / 'script1_outlines' / f'{tag}.geojson',
             'cell_geometry': run / 'cell_geojson' / f'{tag}.geojson'}
    digests = {}
    for role, path in paths.items():
        folder = 'cell_geojson' if role == 'cell_geometry' else path.parent.name
        digests[role] = verify_recorded_source(path, f'{folder}/{path.name}', recorded, sources)
    final_path, reference_path = run / 'islets' / f'{tag}.geojson', reference / 'islets' / f'{tag}.geojson'
    expected_geometry = evidence['geometry'][tag].get('base_file_sha256', '')
    if not re.fullmatch(r'[0-9a-f]{64}', expected_geometry):
        raise ValueError(f'{tag}: missing pinned post-snap geometry hash')
    for path in (final_path, reference_path):
        sources[str(path)] = sha256(path)
        if sources[str(path)] != expected_geometry:
            raise ValueError(f'{path}: final object geometry differs from the preserved generation')
    final = [f for f in features(final_path) if f.get('properties', {}).get('classification', {}).get('name') in CLASSES]
    originals = [f for f in features(paths['original_outlines']) if f.get('properties', {}).get('classification', {}).get('name') == 'Islet']
    if len(final) != evidence['geometry'][tag].get('endocrine_objects'):
        raise ValueError(f'{tag}: pinned final object inventory count differs')
    original_ids = [f.get('id') for f in originals]
    if not originals or any(not x for x in original_ids) or len(set(original_ids)) != len(original_ids):
        raise ValueError(f'{tag}: missing or duplicate original object UUIDs')
    original_index = {key: i for i, key in enumerate(original_ids)}
    original_geometries = [shape(f['geometry']) for f in originals]
    final_geometries = [shape(f['geometry']) for f in final]
    if any(g.is_empty or not g.is_valid or g.area <= 0 for g in original_geometries + final_geometries):
        raise ValueError(f'{tag}: invalid object geometry; no automatic repair is permitted')
    parameters = json.loads(paths['fit'].read_text())
    if parameters.get('tag') != tag or parameters.get('resolution_um') != .7548:
        raise ValueError(f'{tag}: unsupported original classifier grid')
    thresholds = parameters.get('thresholds', {})
    if (set(thresholds) != set(MARKERS) or not any(thresholds.get(m) is not None for m in MARKERS[:3])
            or any(value is not None and (not np.isfinite(value) or value <= 0) for value in thresholds.values())):
        raise ValueError(f'{tag}: incomplete original classifier marker parameters')
    mask = np.load(paths['mask'], mmap_mode='r', allow_pickle=False)
    if mask.ndim != 2 or mask.dtype.kind not in 'bu':
        raise ValueError(f'{tag}: invalid preserved classifier mask')
    support = [original_mask_support(f['geometry'], mask) for f in originals]
    tree = shapely.STRtree(original_geometries)
    single_bounds, single_matches = {}, {}
    for f, g in zip(final, final_geometries):
        if f['properties']['classification']['name'] == 'EndocrineSingle':
            single_bounds.setdefault(tuple(g.bounds), []).append((f['id'], g))
            single_matches[f['id']] = []
    seen_cells = set()
    for cell in streamed_cell_features(paths['cell_geometry']):
        key = cell.get('id')
        if not key or key in seen_cells:
            raise ValueError(f'{tag}: missing or duplicate source cell geometry UUID')
        seen_cells.add(key)
        candidates = single_bounds.get(geometry_bounds(cell['geometry']), [])
        if not candidates:
            continue
        g = shape(cell['geometry'])
        if not g.is_valid or not cell.get('nucleusGeometry') or cell.get('properties', {}).get('objectType') != 'cell':
            continue
        for single_id, single_geometry in candidates:
            if g.equals(single_geometry):
                single_matches[single_id].append(key)
    if seen_cells != set(rs.object_id):
        raise ValueError(f'{tag}: cell geometry and REDSEA UUID inventories differ')
    eligible = rs.set_index('object_id').qc_analysis_eligible
    memberships = ctx.groupby('endocrine_object_id').object_id.agg(set).to_dict()
    cell_context = ctx.set_index('object_id').endocrine_object_id
    matched_owners = {}
    for object_id, matches in single_matches.items():
        for cell_id in matches:
            matched_owners.setdefault(cell_id, []).append(object_id)
    rows, failures = [], []
    for f, g in zip(final, final_geometries):
        object_id = f['id']; cls = f['properties']['classification']['name']
        retained = object_id in original_index
        candidates = [original_index[object_id]] if retained else sorted(tree.query(g, predicate='intersects'))
        overlaps = [(i, float(g.intersection(original_geometries[i]).area)) for i in candidates]
        overlaps = [(i, area) for i, area in overlaps if area > 0]
        supported = [(i, area) for i, area in overlaps if support[i] > 0]
        reasons = []
        if not supported:
            reasons.append('no_positive_area_link_to_supported_original_classifier_object')
        record = {'endocrine_object_id': object_id, 'classifier_lineage_passes': bool(supported),
                  'original_source_relation': 'retained_source_uuid' if retained else 'split_descendant_spatial_support',
                  'original_parent_unique': len(overlaps) == 1,
                  'original_parent_candidates': json.dumps([original_ids[i] for i, _ in overlaps]),
                  'original_supported_parent_candidates': json.dumps([original_ids[i] for i, _ in supported]),
                  'original_overlap_area_pixels': sum(area for _, area in overlaps),
                  'original_max_final_overlap_fraction': max((area / g.area for _, area in overlaps), default=0),
                  'original_classifier_foreground_pixels': sum(support[i] for i, _ in supported),
                  'single_associated_cell_id': None, 'single_association_match_count': None,
                  'single_associated_analysis_eligible': None, 'single_association_matches_centroid_membership': None}
        if cls == 'EndocrineSingle':
            matches = single_matches[object_id]
            if f['properties'].get('measurements', {}).get('Whole cells') != 1:
                reasons.append('single_original_whole_cells_not_one')
            if len(matches) != 1:
                reasons.append('single_does_not_equal_exactly_one_nucleated_cell_outline')
            elif len(matched_owners[matches[0]]) != 1:
                reasons.append('single_source_cell_has_multiple_single_owners')
            record['single_association_match_count'] = len(matches)
            if len(matches) == 1:
                cell_id = matches[0]
                record.update(single_associated_cell_id=cell_id,
                              single_associated_analysis_eligible=bool(eligible.loc[cell_id]),
                              single_associated_context_object_id=cell_context.loc[cell_id],
                              single_association_matches_centroid_membership=memberships.get(object_id, set()) == {cell_id})
        record['classifier_lineage_passes'] = not reasons
        record['classifier_lineage_failures'] = ';'.join(reasons)
        if reasons:
            failures.append({'endocrine_object_id': object_id, 'reasons': reasons})
        rows.append(record)
    frame = pd.DataFrame(rows)
    singles = frame.single_association_match_count.notna()
    summary = {'image': tag, 'final_islets_sha256': expected_geometry,
               'reference_islets_sha256': sources[str(reference_path)],
               'cell_geometry_sha256': digests['cell_geometry'], 'original_outlines_sha256': digests['original_outlines'],
               'fit_sha256': digests['fit'], 'mask_sha256': digests['mask'], 'objects': len(frame),
               'objects_with_classifier_support': int(frame.original_classifier_foreground_pixels.gt(0).sum()),
               'singles': int(singles.sum()),
               'singles_with_exact_cell_association': int((singles & frame.single_association_match_count.eq(1)).sum()),
               'singles_with_centroid_membership_difference': int((singles & frame.single_association_matches_centroid_membership.eq(False)).sum()),
               'split_descendants': int(frame.original_source_relation.eq('split_descendant_spatial_support').sum()),
               'split_descendants_with_multiple_parent_candidates': int((frame.original_source_relation.eq('split_descendant_spatial_support') & ~frame.original_parent_unique).sum()),
               'failures': failures}
    return frame, summary


def summarize(frame: pd.DataFrame) -> dict:
    negative = frame.all_analysis_cells_negative
    return {'objects': len(frame),
            'objects_with_analysis_cells': int(frame.analysis_cells.gt(0).sum()),
            'objects_without_context_cells': int(frame.context_cells.eq(0).sum()),
            'objects_without_analysis_cells': int(frame.no_analysis_cells.sum()),
            'single_membership_anomalies': int(frame.single_membership_anomaly.sum()),
            'membership_anomalies': int(frame.membership_anomaly.sum()),
            'objects_with_any_endocrine_positive_cell': int(frame.any_endocrine_positive_cells.gt(0).sum()),
            'objects_with_all_analysis_cells_negative': int(negative.sum()),
            'context_cells': int(frame.context_cells.sum()),
            'analysis_cells': int(frame.analysis_cells.sum()),
            'any_endocrine_positive_cells': int(frame.any_endocrine_positive_cells.sum()),
            'all_endocrine_negative_cells': int((frame.analysis_cells - frame.any_endocrine_positive_cells).sum()),
            'analysis_cells_in_fully_negative_objects': int(frame.loc[negative, 'analysis_cells'].sum())}


def audit(run: Path, out: Path, *, levels_path: Path | None = None,
          images: list[str] | None = None, classifier_run: Path | None = None,
          classifier_replay: Path | None = None, classifier_replay_sha256: str | None = None,
          geometry_reference: Path | None = None, geometry_provenance: Path | None = None,
          geometry_provenance_sha256: str | None = None) -> dict:
    """Write a read-only audit and return its report, including ``passes``.

    ``levels_path`` lets levels.py audit its new candidate table, including a
    --report trial, without accidentally reading the run's previous levels.
    """
    run, out = run.resolve(), out.resolve()
    if out.exists():
        raise FileExistsError(f'{out}: output already exists; choose a fresh directory')
    levels_path = (levels_path or run / 'typing/levels.csv').resolve()
    tags = images if images is not None else sorted(p.stem for p in (run / 'islets').glob('*.geojson'))
    if not tags or len(tags) != len(set(tags)):
        raise ValueError('images must be nonempty and unique')
    lineage_arguments = (classifier_run, classifier_replay, classifier_replay_sha256,
                         geometry_reference, geometry_provenance, geometry_provenance_sha256)
    if any(x is not None for x in lineage_arguments) and not all(x is not None for x in lineage_arguments):
        raise ValueError('Classifier-lineage mode requires all six explicit source paths and SHA256 pins')
    lineage = (load_lineage_evidence(*lineage_arguments, tags)
               if all(x is not None for x in lineage_arguments) else None)
    levels = pd.read_csv(levels_path)
    levels = levels.loc[levels.image.isin(tags) & levels.marker.isin(MARKERS)]
    if levels.duplicated(['image', 'marker']).any():
        raise ValueError(f'{levels_path}: duplicate endocrine levels')
    sources = {str(levels_path): sha256(levels_path)}
    if lineage:
        for key in ('classifier_replay', 'geometry_provenance'):
            sources[lineage[key]['path']] = lineage[key]['sha256']
    frames, checks, lineage_checks = [], [], []
    for tag in tags:
        selected = levels.loc[levels.image.eq(tag)].set_index('marker')
        if set(selected.index) != set(MARKERS):
            raise ValueError(f'{tag}: all four endocrine marker levels are required')
        if not selected.measurement.eq('whole_cell').all():
            raise ValueError(f'{tag}: endocrine positivity must use whole-cell measurements')
        if not np.isfinite(selected.level_log.to_numpy(float)).all():
            raise ValueError(f'{tag}: nonfinite endocrine positivity level')
        rs_path = run / 'typing/redsea' / f'{tag}.parquet'
        ctx_path = run / 'typing/context' / f'{tag}.parquet'
        obj_path = run / 'islets' / f'{tag}.geojson'
        for source in (rs_path, ctx_path, obj_path):
            sources[str(source)] = sha256(source)
        obj = object_inventory(obj_path)
        rs = pd.read_parquet(rs_path, columns=['object_id', 'qc_analysis_eligible'] + ['Cell__' + m for m in MARKERS])
        ctx = pd.read_parquet(ctx_path, columns=['object_id', 'qc_analysis_eligible',
                                               'endocrine_object_id', 'endocrine_object_class'])
        for name, frame in (('REDSEA', rs), ('context', ctx)):
            if frame.object_id.isna().any() or frame.object_id.duplicated().any():
                raise ValueError(f'{tag}: missing or duplicate {name} cell IDs')
            if not frame.qc_analysis_eligible.isin([False, True]).all():
                raise ValueError(f'{tag}: invalid {name} analysis-eligibility flags')
        d = rs.merge(ctx, on='object_id', how='outer', validate='one_to_one',
                     suffixes=('_redsea', '_context'), indicator=True)
        if not d._merge.eq('both').all():
            raise ValueError(f'{tag}: REDSEA and context cell universes differ')
        if not d.qc_analysis_eligible_redsea.eq(d.qc_analysis_eligible_context).all():
            raise ValueError(f'{tag}: REDSEA and context analysis eligibility disagree')
        in_object = d.endocrine_object_id.notna()
        if not in_object.eq(d.endocrine_object_class.notna()).all():
            raise ValueError(f'{tag}: incomplete context object membership')
        class_map = obj.set_index('endocrine_object_id').endocrine_object_class
        if not d.loc[in_object, 'endocrine_object_class'].eq(
                d.loc[in_object, 'endocrine_object_id'].map(class_map)).all():
            raise ValueError(f'{tag}: context references an absent object or mismatched object class')
        d = d.loc[in_object & d.qc_analysis_eligible_redsea.astype(bool)].copy()
        for marker in MARKERS:
            x = d['Cell__' + marker].to_numpy(float)
            d['pos__' + marker] = np.log1p(np.clip(np.nan_to_num(x, nan=0.0), 0, None)) > selected.loc[marker, 'level_log']
        d['any_endocrine_positive'] = d[['pos__' + m for m in MARKERS]].any(axis=1)
        count_columns = {m + '_positive_cells': ('pos__' + m, 'sum') for m in MARKERS}
        count_columns.update({'analysis_cells': ('object_id', 'size'),
                              'any_endocrine_positive_cells': ('any_endocrine_positive', 'sum')})
        counts = d.groupby('endocrine_object_id').agg(**count_columns)
        context_counts = ctx.groupby('endocrine_object_id').size().rename('context_cells')
        result = obj.merge(counts, on='endocrine_object_id', how='left', validate='one_to_one')
        result = result.merge(context_counts, on='endocrine_object_id', how='left', validate='one_to_one')
        for column in list(count_columns) + ['context_cells']:
            result[column] = result[column].fillna(0).astype(int)
        result['all_analysis_cells_negative'] = result.analysis_cells.gt(0) & result.any_endocrine_positive_cells.eq(0)
        result['no_analysis_cells'] = result.analysis_cells.eq(0)
        result['single_membership_anomaly'] = result.endocrine_object_class.eq('EndocrineSingle') & result.analysis_cells.ne(1)
        result['membership_anomaly'] = result.no_analysis_cells | result.single_membership_anomaly
        result['image'] = tag
        if lineage:
            proof, lineage_check = audit_classifier_lineage(run, tag, lineage, rs, ctx, sources)
            result = result.merge(proof, on='endocrine_object_id', how='left', validate='one_to_one')
            if result.classifier_lineage_passes.isna().any():
                raise ValueError(f'{tag}: incomplete classifier-lineage object inventory')
            lineage_checks.append(lineage_check)
        frames.append(result)
        checks.append({'image': tag, 'redsea_cells': len(rs), 'context_cells': len(ctx),
                       'analysis_cells': int(rs.qc_analysis_eligible.sum()), 'analysis_cells_in_objects': len(d),
                       'exact_context_join': True, 'matching_analysis_eligibility': True,
                       'object_class_matches_geometry': True})
    all_objects = pd.concat(frames, ignore_index=True)
    totals = summarize(all_objects)
    corrected_passes = not bool(all_objects.all_analysis_cells_negative.any() or all_objects.membership_anomaly.any())
    report = {'schema': SCHEMA, 'run': str(run), 'levels': str(levels_path), 'images': tags,
              'audit_mode': 'classifier_lineage' if lineage else 'corrected_cell',
              'passes': bool(all_objects.classifier_lineage_passes.all()) if lineage else corrected_passes,
              'corrected_cell_diagnostic_passes': corrected_passes,
              'rule': 'Every final endocrine object has at least one analysis-eligible cell positive for INS, GCG, SST or CHGA; '
                      'EndocrineSingle has exactly one analysis-eligible member.',
              'totals': totals,
              'by_class': {label: summarize(group) for label, group in all_objects.groupby('endocrine_object_class')},
              'by_image_and_class': [dict(image=tag, object_class=label, **summarize(group))
                                     for (tag, label), group in all_objects.groupby(['image', 'endocrine_object_class'])],
              'checks': checks, 'source_sha256': sources,
              'limitations': ['This checks object/cell consistency; it does not replay the raw-image pixel classifier.',
                              'Other cells within a marker-positive object need not be endocrine-marker positive.',
                              'Context membership is centroid-based and selects the smallest containing object.']}
    if lineage:
        report['corrected_cell_diagnostic_rule'] = report['rule']
        report['rule'] = ('Every preserved final endocrine object is linked to a supported original classifier region; '
                          'every Single equals exactly one nucleated source cell outline. Centroid membership '
                          'and corrected-cell positivity are independent retained diagnostics.')
        count_fields = ('objects', 'objects_with_classifier_support', 'singles',
                        'singles_with_exact_cell_association', 'singles_with_centroid_membership_difference',
                        'split_descendants', 'split_descendants_with_multiple_parent_candidates')
        report['classifier_lineage'] = {
            'passes': report['passes'], 'classifier_replay': lineage['classifier_replay'],
            'geometry_provenance': lineage['geometry_provenance'], 'images': lineage_checks,
            'totals': {**{key: sum(row[key] for row in lineage_checks) for key in count_fields},
                       'failures': sum(len(row['failures']) for row in lineage_checks)}}
        report['limitations'] = [
            'Classifier support validates original region provenance, not a hormone identity for each cell.',
            'A negative corrected whole-cell mean is not evidence that the original classifier region lacks signal.',
            'Single association is an independent exact outline relationship; centroid memberships are not rewritten.',
            'Split descendants retain all supported spatial parent candidates; no unique historical parent is invented.',
            'No minimum source-overlap fraction is added; small preserved overlaps remain reported for review.',
            'Passing this object audit alone does not establish cell-typing accuracy or authorize deployment.']
    out.mkdir(parents=True)
    all_objects.to_csv(out / 'objects.csv', index=False)
    all_objects.loc[all_objects.all_analysis_cells_negative].to_csv(out / 'fully_negative_objects.csv', index=False)
    all_objects.loc[all_objects.membership_anomaly].to_csv(out / 'membership_anomalies.csv', index=False)
    if lineage:
        all_objects.loc[~all_objects.classifier_lineage_passes].to_csv(out / 'classifier_lineage_failures.csv', index=False)
        all_objects.loc[all_objects.endocrine_object_class.eq('EndocrineSingle')].to_csv(
            out / 'single_associations.csv', index=False)
    (out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--levels', type=Path, help='candidate levels table; default: RUN/typing/levels.csv')
    parser.add_argument('--images', nargs='+', help='default: every image in the final islets directory')
    parser.add_argument('--classifier-run', type=Path, help='original classifier inputs; activates explicit lineage mode')
    parser.add_argument('--classifier-replay', type=Path, help='all-20 exact classifier replay receipt')
    parser.add_argument('--classifier-replay-sha256', help='required explicit SHA256 pin of the replay receipt')
    parser.add_argument('--geometry-reference', type=Path, help='preserved original post-snap run')
    parser.add_argument('--geometry-provenance', type=Path, help='all-20 historical final-geometry hash receipt')
    parser.add_argument('--geometry-provenance-sha256', help='required explicit SHA256 pin of geometry provenance')
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output already exists; choose a fresh directory')
    report = audit(args.run, args.out, levels_path=args.levels, images=args.images,
                   classifier_run=args.classifier_run, classifier_replay=args.classifier_replay,
                   classifier_replay_sha256=args.classifier_replay_sha256,
                   geometry_reference=args.geometry_reference, geometry_provenance=args.geometry_provenance,
                   geometry_provenance_sha256=args.geometry_provenance_sha256)
    print(json.dumps({'passes': report['passes'], 'totals': report['totals'], 'by_class': report['by_class']}, indent=2))
    return 0 if report['passes'] else 4


if __name__ == '__main__':
    sys.exit(main())
