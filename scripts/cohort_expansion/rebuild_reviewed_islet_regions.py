#!/usr/bin/env python3
"""Recompute regional measurements only for a hash-bound, owner-approved donor.

The frozen cohort and raw cells remain inputs. Geometry, membership, and regional
outputs are new files; boundary approval never certifies phenotype gates.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'research'))
from build_candidate import process
from common import MARKERS, read_cohort, sha256, write_json
from phenotypes import read_rules


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('input-root', 'cohort', 'project-root', 'export-root', 'exact-geometry-root', 'gate-registry', 'output-root'):
        parser.add_argument('--' + arg, type=Path, required=True)
    args = parser.parse_args()
    root = args.input_root.resolve(strict=True)
    project = args.project_root.resolve(strict=True)
    exports = args.export_root.resolve(strict=True)
    out = args.output_root.resolve()
    approval_file = project / 'boundary_approval.json'
    approval = json.loads(approval_file.read_text())
    promotion_file = exports / 'promotion_manifest.json'
    promotion = json.loads(promotion_file.read_text())
    cid = str(approval['case_id'])
    if (approval['decision'] != 'owner_visually_approved'
            or promotion['status'] != 'approved_boundaries_promoted_and_native_roundtrip_passed'
            or promotion['case_id'] != cid or promotion['approval_sha256'] != sha256(approval_file)):
        raise ValueError('An explicit, verified boundary promotion is required')
    geometry = exports / 'annotations.geojson.gz'
    accepted = exports / 'accepted_islets.geojson'
    if (sha256(geometry) != promotion['geometry_sha256']
            or sha256(accepted) != promotion['core_geometry_sha256']):
        raise ValueError('Approved geometry exports changed')
    staging = json.loads((project / 'staging_manifest.json').read_text())
    native = Path(next(d['review_data'] for d in staging['donors'] if d['case_id'] == cid))
    if sha256(native) != promotion['output_native_sha256']:
        raise ValueError('Promoted native data changed')
    exact_root = args.exact_geometry_root.resolve(strict=True)
    exact_file = exact_root / 'verification.json'
    exact = json.loads(exact_file.read_text())
    geometry = exact_root / 'annotations.geojson.gz'
    if (exact['status'] != 'exact_native_polygon_vertices_exported'
            or exact['native_data_sha256'] != promotion['output_native_sha256']
            or exact['geometry_sha256'] != sha256(geometry)
            or exact['cores'] != promotion['accepted_islets']):
        raise ValueError('Exact native polygon geometry provenance differs')
    cohort = read_cohort(args.cohort, root)
    selected = cohort.loc[cohort.case_id == cid].copy()
    if len(selected) != 1:
        raise ValueError('Approved donor missing from explicit cohort')
    original_geometry = selected.geometry.iloc[0]
    selected['geometry'] = str(geometry)
    selected['annotation_review'] = 'owner_visually_approved_20260918'
    registry = pd.read_csv(args.gate_registry, dtype={'case_id': str})
    if registry.duplicated(['case_id', 'marker']).any() or not registry.validated.isin([True, False]).all():
        raise ValueError('Gate provenance must be unique and explicitly Boolean')
    gates = registry.loc[registry.case_id == cid].rename(
        columns={'case_id': 'imageid', 'marker': 'markers', 'gate_log1p': 'gate'})
    if gates.empty:
        raise ValueError('Explicit gate provenance is required, including unvalidated gates')
    out.mkdir(parents=True, exist_ok=False)
    for sub in ('analysis', 'raw'):
        (out / sub).mkdir()
    selected.to_csv(out / 'source_manifest.csv', index=False)
    gates.to_csv(out / 'gate_provenance.csv', index=False)
    identity = dict(case_id=cid, original_geometry=original_geometry,
        approved_geometry=str(geometry), approved_geometry_sha256=sha256(geometry),
        raw_cells=selected.cells.iloc[0], raw_cells_sha256=sha256(selected.cells.iloc[0]),
        cohort_sha256=sha256(args.cohort), gate_registry_sha256=sha256(args.gate_registry),
        approval_sha256=sha256(approval_file), promotion_manifest_sha256=sha256(promotion_file),
        exact_geometry_verification_sha256=sha256(exact_file),
        runner_sha256=sha256(__file__),
        regional_builder_sha256=sha256(Path(__file__).parents[1] / 'research/build_candidate.py'),
        phenotype_rules_sha256=sha256(root / 'data/phenotype_rules.csv'),
        phenotype_engine_sha256=sha256(Path(__file__).parents[1] / 'research/phenotypes.py'),
        boundary_approval_does_not_validate_phenotypes=True)
    write_json(out / 'build_identity.json', identity)
    core_types = exact['core_roi_type_by_name']
    if len(core_types) != promotion['accepted_islets'] or set(core_types.values()) - {'PolygonROI', 'GeometryROI'}:
        raise ValueError('Unsupported or incomplete native containment types')
    process(next(selected.itertuples(index=False)), out,
            read_rules(root / 'data/phenotype_rules.csv'), gates, gates, core_roi_types=core_types)
    membership = pd.read_parquet(out / 'analysis' / f'{cid}-memberships.parquet')
    native_members = pd.read_csv(exports / 'native_core_memberships.csv', dtype=str)
    native_members['islet_id'] = cid + '_' + native_members.islet_key
    expected = set(zip(native_members.islet_id, native_members.cell_id))
    actual = set(zip(membership.islet_id, membership.object_id))
    if actual != expected or membership.object_id.duplicated().any():
        write_json(out / 'membership_discrepancy.json', dict(
            only_native=sorted(expected - actual), only_offline=sorted(actual - expected)))
        raise ValueError('Offline core membership does not match native QuPath containment')
    regional = pd.read_parquet(out / 'analysis' / f'{cid}.parquet')
    core = regional.query("compartment == 'core' and band_um == 20").set_index('islet_key')
    features = json.loads(accepted.read_text())['features']
    for feature in features:
        p = feature['properties']
        row, measures = core.loc[p['name']], p['measurements']
        for key, native_key in [('raw_cells', 'raw cells'), ('measured_cells', 'measured cells')]:
            if row[key] != measures['Rebuilt core: ' + native_key]:
                raise ValueError('Native/offline core counts differ: ' + p['name'])
        a = np.asarray([row['median_log1p_' + m] for m in MARKERS], float)
        b = np.asarray([measures.get('Rebuilt core: median log1p ' + m, np.nan) for m in MARKERS], float)
        if not np.allclose(a, b, equal_nan=True, atol=2e-6, rtol=0):
            raise ValueError('Native/offline marker medians differ: ' + p['name'])
    # Source cell exports retain historical parents. Keep them as archival
    # fields and join the newly verified core assignment explicitly.
    cell_file = out / 'raw' / f'{cid}-cells.parquet'
    cells = pd.read_parquet(cell_file).rename(columns={'Parent': 'archival_parent'})
    cells = cells.merge(membership, on='object_id', how='left', validate='one_to_one')
    cells = cells.rename(columns={'islet_id': 'current_islet_id'})
    cells.to_parquet(cell_file, index=False)
    regional.to_csv(out / 'analysis' / f'{cid}-regions.csv.gz', index=False)
    core.reset_index().to_csv(out / 'analysis' / f'{cid}-core_summary.csv', index=False)
    audit = pd.read_csv(out / 'analysis' / f'{cid}-phenotype-audit.csv')
    complete = json.loads((out / 'analysis' / f'{cid}-complete.json').read_text())
    result = dict(status='approved_geometry_regions_verified', case_id=cid,
        **{k: complete[k] for k in ('raw_cells', 'measured_cells', 'unmeasured_cells', 'islets')},
        native_core_cells=len(native_members), exact_native_offline_membership_match=True,
        core_coordinate_source='Native ROI vertices without export rounding, with type-specific boundary containment',
        all_31_native_offline_core_marker_summaries_match=True,
        compartments=['core', 'peri'], peri_bands_um=[10, 20, 30],
        primary_spatial_assessable_islets=int((core.geometry_valid & ~core.overlap).sum()),
        overlapping_30um_footprints=int(core.overlap.sum()),
        minimum_10_measured_core_cells=int((core.measured_cells >= 10).sum()),
        validated_phenotypes=audit.loc[audit.validated, 'phenotype'].tolist(),
        unavailable_phenotypes=audit.loc[~audit.validated, 'phenotype'].tolist(),
        structural_density='union of positive mask area / valid compartment area',
        native_expansions='core plus tissue-clipped 20 um; offline peri regions subtract all cores',
        acceptance='HPAP boundaries accepted; phenotype gate approval remains separate',
        artifacts={str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob('*')) if p.is_file()})
    write_json(out / 'verification.json', result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
