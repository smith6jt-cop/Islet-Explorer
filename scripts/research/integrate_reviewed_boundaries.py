#!/usr/bin/env python3
"""Build an independent app revision from verified, owner-approved HPAP geometry.

The other twenty donors keep their explicit existing inputs and review status.
Never project replacement islets onto an axis fitted from rejected geometry.
"""
import argparse
import gzip
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from common import MARKERS, ENDOCRINE, sha256, write_json, require_cohort
from trajectory import apply_identity_review


def replace_donor(old, new, column='case_id'):
    return pd.concat([old.loc[old[column].astype(str) != '78'], new], ignore_index=True)


def features_from_regions(regions):
    core = regions.query("band_um == 20 and compartment == 'core'").copy()
    peri = regions.query("band_um == 20 and compartment == 'peri'").copy()
    columns = ['median_log1p_' + m for m in MARKERS]
    features = core.drop(columns=columns).merge(core[['islet_id'] + columns].rename(
        columns={c: 'core__' + c[13:] for c in columns}), on='islet_id', validate='one_to_one')
    features = features.merge(peri[['islet_id'] + columns].rename(
        columns={c: 'peri__' + c[13:] for c in columns}), on='islet_id', validate='one_to_one')
    features['eligible'] = features.geometry_valid & (features.measured_cells >= 10) & (features.core_tissue_coverage > .99)
    features['exclusion_reason'] = np.where(~features.geometry_valid, 'invalid_geometry',
        np.where(features.measured_cells < 10, 'fewer_than_10_measured_core_cells',
                 np.where(features.core_tissue_coverage <= .99, 'core_outside_tissue', '')))
    missing = ~np.isfinite(features[[c for c in features if c.startswith(('core__', 'peri__'))]]).all(axis=1)
    features.loc[missing, 'eligible'] = False
    features.loc[missing, 'exclusion_reason'] = 'empty_or_unmeasured_compartment'
    return features


def runtime_summaries(source, out, regions):
    """Replace every old HPAP row; unreviewed positive counts remain missing."""
    frames = {p.stem: pd.read_parquet(p) for p in (source / 'summary').glob('*.parquet')}
    core = regions.query("band_um == 20 and compartment == 'core'").set_index('islet_key')
    peri = regions.query("band_um == 20 and compartment == 'peri'").set_index('islet_key')
    identity = frames['islets'].loc[frames['islets']['Case ID'].astype(str).eq('78')].iloc[0]
    varying = {'islet_key', 'cells_total', 'Ins_single', 'Glu_single', 'Stt_single', 'Multi_Pos',
               'Triple_Neg', 'Ins_any', 'Glu_any', 'Stt_any', 'islet_diam_um', 'core_region_um2', 'peri_region_um2'}
    metadata = {k: v for k, v in identity.items() if k not in varying}
    area_rows, comp_rows, target_rows, marker_rows = [], [], [], []
    for name, c in core.iterrows():
        p = peri.loc[name]
        base = dict(metadata, islet_key=name, islet_diam_um=float(np.sqrt(4*c.valid_area_um2/np.pi)))
        area_rows.append({k: base[k] for k in frames['core_area'].columns})
        comp = dict(base, cells_total=int(c.raw_cells), core_region_um2=c.valid_area_um2, peri_region_um2=p.valid_area_um2)
        for key in varying - {'islet_key', 'cells_total', 'islet_diam_um', 'core_region_um2', 'peri_region_um2'}:
            comp[key] = np.nan
        comp_rows.append(comp)
        for kind, parts in [('islet_core', [c]), ('islet_band', [p]), ('islet_union', [c, p])]:
            region_area = sum(r.valid_area_um2 for r in parts)
            for structure in ['Nerve', 'Capillary', 'Lymphatic']:
                positive = sum(r[structure.lower() + '_area_um2'] for r in parts)
                target_rows.append(dict(base, type=kind, **{'class': structure}, area_um2=positive,
                    region_um2=region_area, area_density=positive/region_area if region_area > 0 else np.nan, count=np.nan))
            for marker in MARKERS:
                marker_rows.append(dict(base, region_type=kind, marker=marker,
                    n_cells=sum(r.measured_cells for r in parts), pos_count=np.nan, pos_frac=np.nan, pos_pct=np.nan))
    new = dict(core_area=area_rows, islets=comp_rows, targets=target_rows, markers=marker_rows)
    (out / 'summary').mkdir()
    for name, original in frames.items():
        replacement = pd.DataFrame(new[name]).reindex(columns=original.columns) if name in new else original.iloc[:0]
        frame = replace_donor(original, replacement, 'Case ID')
        frame.to_parquet(out / 'summary' / (name + '.parquet'), index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['input-root', 'output-root', 'review-project', 'regional-root', 'exact-geometry-root']:
        parser.add_argument('--' + key, required=True, type=Path)
    args = parser.parse_args()
    source, out = args.input_root.resolve(strict=True), args.output_root.resolve()
    project = args.review_project.resolve(strict=True)
    regional = args.regional_root.resolve(strict=True)
    exact_root = args.exact_geometry_root.resolve(strict=True)
    delivery = json.loads((project / 'delivery_manifest.json').read_text())
    verification = json.loads((regional / 'verification.json').read_text())
    exact = json.loads((exact_root / 'verification.json').read_text())
    geometry = exact_root / 'annotations.geojson.gz'
    if delivery['accepted_boundaries'] != 1431 or verification['status'] != 'approved_geometry_regions_verified':
        raise ValueError('Verified owner-approved HPAP delivery required')
    if sha256(geometry) != exact['geometry_sha256']:
        raise ValueError('Verified geometry changed')
    for rel, expected in verification['artifacts'].items():
        if sha256(regional / rel) != expected:
            raise ValueError('Verified regional artifact changed: ' + rel)
    cohort = pd.read_csv(source / 'cohort.csv', dtype=str)
    if len(cohort) != 21 or set(cohort.case_id) != set(pd.read_csv(source / 'donor_metadata.csv', dtype=str).case_id):
        raise ValueError('Explicit 21-donor source required')
    out.mkdir(parents=True, exist_ok=False)
    for folder in ['images', 'raw', 'gson', 'regional', 'analysis', 'geometry', 'cells']:
        for path in sorted((source / folder).rglob('*')):
            if not path.is_file():
                continue
            if folder == 'cells' and path.name.startswith('78_'):
                continue
            if folder in {'gson', 'geometry'} and path.name.startswith('0078.'):
                continue
            if folder == 'regional' and '_0078.' in path.name:
                continue
            if folder == 'raw' and path.name == '78-cells.parquet':
                continue
            if folder == 'analysis' and (path.name.startswith('78.') or path.name.startswith('78-') or path.name in {'regions.parquet', 'phenotype_audit.csv'}):
                continue
            dest = out / path.relative_to(source)
            dest.parent.mkdir(parents=True, exist_ok=True)
            os.link(path, dest)
    for folder in ['logs', 'evidence', 'trajectory', 'provenance', 'gson', 'cells']:
        (out / folder).mkdir(exist_ok=True)
    for name in ['donor_metadata.csv', 'phenotype_rules.csv']:
        shutil.copy2(source / name, out / name)
    for name in ['cohort.csv', 'source_manifest.csv']:
        frame = pd.read_csv(source / name, dtype=str)
        frame.loc[frame.case_id.eq('78'), 'geometry'] = str(geometry)
        frame.loc[frame.case_id.eq('78'), 'regional_root'] = str(regional / 'analysis')
        frame.loc[frame.case_id.eq('78'), 'annotation_review'] = 'owner_visually_approved_20260918'
        frame.to_csv(out / name, index=False)
    shutil.copy2(geometry, out / 'gson/0078.geojson.gz')
    shutil.copy2(exact_root / 'verification.json', out / 'provenance/78-native-geometry.json')
    shutil.copy2(project / 'boundary_approval.json', out / 'provenance/78-boundary-approval.json')
    shutil.copy2(regional / 'raw/78-cells.parquet', out / 'raw/78-cells.parquet')
    for path in (regional / 'analysis').iterdir():
        shutil.copy2(path, out / 'analysis' / path.name)
    regions = pd.read_parquet(regional / 'analysis/78.parquet')
    all_regions = replace_donor(pd.read_parquet(source / 'analysis/regions.parquet'), regions)
    all_regions.to_parquet(out / 'analysis/regions.parquet', index=False)
    new_features = features_from_regions(regions)
    features = replace_donor(pd.read_csv(source / 'trajectory/features.csv', dtype={'case_id': str}), new_features)
    features['boundary_review_status'] = np.where(features.case_id.eq('78'), 'owner_approved', 'review_pending')
    features.to_csv(out / 'trajectory/features.csv', index=False)
    counts = replace_donor(pd.read_csv(source / 'trajectory/counts.csv', dtype={'case_id': str}),
                           pd.read_parquet(regional / 'analysis/78-counts.parquet'))
    counts.to_csv(out / 'trajectory/counts.csv', index=False)
    audit = replace_donor(pd.read_csv(source / 'trajectory/phenotype_availability.csv', dtype={'case_id': str}),
                          pd.read_csv(regional / 'analysis/78-phenotype-audit.csv', dtype={'case_id': str}))
    audit.to_csv(out / 'trajectory/phenotype_availability.csv', index=False)
    audit.to_csv(out / 'analysis/phenotype_audit.csv', index=False)
    denom = pd.read_csv(source / 'trajectory/endocrine_counts.csv', dtype={'case_id': str})
    hpap_denom = new_features[['islet_id', 'case_id', 'endocrine_count', 'measured_cells', 'geometry_valid']].rename(columns={'measured_cells': 'measured_core_cells'})
    hpap_denom['validated'] = hpap_denom.geometry_valid & hpap_denom.endocrine_count.notna()
    missing = sorted(set(';'.join(audit.loc[audit.case_id.eq('78') & audit.phenotype.isin(ENDOCRINE), 'missing_validated_gates'].fillna('')).split(';')) - {''})
    hpap_denom['reason'] = 'gate_review_required: ' + ','.join(missing)
    hpap_denom['definition'] = 'unique core alpha/beta/delta cell union'
    hpap_denom['rule_sha256'] = sha256(out / 'phenotype_rules.csv')
    replace_donor(denom, hpap_denom).to_csv(out / 'trajectory/endocrine_counts.csv', index=False)
    shutil.copy2(source / 'trajectory/marker_panel.json', out / 'trajectory/marker_panel.json')
    cm = json.loads((source / 'trajectory/count_manifest.json').read_text())
    cm.update(cohort_sha256=sha256(out / 'cohort.csv'), counts_sha256=sha256(out / 'trajectory/counts.csv'),
              endocrine_counts_sha256=sha256(out / 'trajectory/endocrine_counts.csv'), geometry_revision='HPAP owner-approved native outlines')
    write_json(out / 'trajectory/count_manifest.json', cm)
    lookup = features[['case_id', 'islet_key', 'islet_id', 'centroid_x_um', 'centroid_y_um', 'area_um2']].rename(columns={'islet_id': 'combined_islet_id'})
    lookup.to_csv(out / 'islet_spatial_lookup.csv', index=False)
    sources = pd.read_csv(out / 'source_manifest.csv', dtype=str)
    hashes = {r.case_id: sha256(r.geometry) for r in sources.itertuples(index=False)}
    review = features[['islet_id', 'case_id']].copy()
    review['decision'] = np.where(review.case_id.eq('78'), 'confirmed_islet', 'unreviewed')
    review['geometry_sha256'] = review.case_id.map(hashes)
    review['review_note'] = np.where(review.case_id.eq('78'), 'Owner visual approval 2026-09-18; exact native geometry and memberships verified', '')
    review.to_csv(out / 'trajectory/islet_review.csv', index=False)
    try:
        apply_identity_review(features, review, hashes)
    except ValueError as error:
        refit_reason = str(error)
    else:
        raise ValueError('Review contract changed: run a full refit instead of emitting unavailable coordinates')
    coords = features.drop(columns=[c for c in features if c.startswith(('core__', 'peri__'))]).copy()
    for name in ['pseudotime_raw', 'pseudotime_percentile'] + ['DC' + str(i) for i in range(1, 11)]:
        coords[name] = np.nan
    coords['supported'] = False; coords['component'] = -1; coords['landmark'] = False
    coords['coordinate_status'] = 'awaiting_reviewed_refit_after_geometry_change'
    coords.to_csv(out / 'trajectory/coordinates.csv', index=False)
    pd.DataFrame(columns=['from', 'to', 'weight']).to_csv(out / 'trajectory/edges.csv', index=False)
    write_json(out / 'trajectory/specification.json', dict(status='awaiting_reviewed_refit', columns=[s+'__'+m for s in ['core','peri'] for m in MARKERS if m not in ['INS','GCG','SST','CHGA']],
        features_sha256=sha256(out / 'trajectory/features.csv'), cohort_sha256=sha256(out / 'cohort.csv'),
        source_geometry_sha256=hashes, islet_review_sha256=sha256(out / 'trajectory/islet_review.csv'),
        root=None, previous_model_reused=False, refit_precondition=refit_reason))
    write_json(out / 'trajectory/validation.json', dict(status='awaiting_reviewed_refit', coordinates_available=False,
        eligible=int(features.eligible.sum()), supported=0, bootstrap_completed=0, lodo_completed=0,
        interpretation='Corrected geometry and measurements are integrated. All-donor pseudotime requires the remaining islet-identity reviews and a new offline fit.'))
    write_json(out / 'trajectory/quality_review.json', dict(status='awaiting_reviewed_refit', current_fit_accepted=False,
        message='HPAP-078 now uses 1,431 owner-approved endocrine boundaries. The previous ordering was invalidated by this geometry change. Remaining donor boundaries require review before a new all-donor fit.',
        approved_cases=['78'], awaiting_review_cases=sorted(set(cohort.case_id)-{'78'})))
    rows=[]
    new_donors={d['case_id']: d for d in delivery['donors']}
    for cid, group in features.groupby('case_id', sort=False):
        record=new_donors.get(cid, {})
        rows.append(dict(case_id=cid, islets=len(group), approved=cid=='78',
            status='Owner-approved corrected boundaries' if cid=='78' else ('Previous boundaries shown; regenerated candidates await visual review' if cid in new_donors else 'Existing reference boundaries; identity review pending'),
            pending_new_candidates=record.get('unreviewed_candidates',0)))
    pd.DataFrame(rows).to_csv(out / 'trajectory/boundary_review.csv', index=False)
    runtime_summaries(source, out, regions)
    shutil.copy2(source / 'cell_lookup.json', out / 'cell_lookup.json')
    build = json.loads((source / 'build_report.json').read_text())
    replacement = json.loads((regional / 'analysis/78-complete.json').read_text())
    build['donor_reports'] = [replacement if str(r['case_id'])=='78' else r for r in build['donor_reports']]
    build.update(islets=len(features), eligible_islets=int(features.eligible.sum()), boundary_revision='owner_approved_HPap78', trajectory_status='awaiting_reviewed_refit')
    write_json(out / 'build_report.json', build)
    write_json(out / 'candidate.json', dict(schema_version='3', donors=21, expected_images=sorted(Path(p).name for p in cohort.image),
        cohort_sha256=sha256(out / 'cohort.csv'), scientific_status='partial_boundary_review', ufrc_accepted=False, research_results_in_app=False))
    write_json(out / 'revision.json', dict(kind='approved_boundary_integration', source_root=str(source),
        source_manifest_sha256=sha256(source / 'manifest.json'), approved_cases=['78'],
        new_geometry_sha256=hashes['78'], regional_verification_sha256=sha256(regional / 'verification.json'),
        previous_trajectory_reused=False, previous_hpap_ids_removed=True,
        instruction='Linked inputs are immutable; replace files atomically, never edit them in place.'))
    print(json.dumps(dict(islets=len(features), hpap_islets=len(new_features), eligible=int(features.eligible.sum()),
                         trajectory='awaiting_reviewed_refit', old_hpap_files_omitted=True)), flush=True)


if __name__ == '__main__':
    main()
