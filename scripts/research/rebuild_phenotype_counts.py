#!/usr/bin/env python3
"""Versioned gate evidence and retained-rule counts, independent of archival labels."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from common import MARKERS, RETAINED, ENDOCRINE, require_cohort, sha256, write_json
from phenotypes import read_rules, audit
from build_candidate import geometry_data, membership


def gate_registry(root, old_gates):
    rows = []
    cohort = pd.read_csv(root/'cohort.csv', dtype=str)
    legacy = pd.read_csv(old_gates, dtype={'imageid': str})
    legacy_hash = sha256(old_gates)
    for cid in cohort.case_id:
        composition = root/'regional'/f'islet_composition_{int(cid):04d}.tsv'
        comp = pd.read_csv(composition, sep='\t')
        thresholds = comp[['Ins_thr', 'Glu_thr', 'Stt_thr']].drop_duplicates()
        if len(thresholds) != 1: raise ValueError(f'Nonconstant exported hormone thresholds: {cid}')
        hormone = dict(zip(['INS', 'GCG', 'SST'], thresholds.iloc[0].to_numpy(float)))
        previous = legacy.loc[legacy.imageid == cid].set_index('markers')
        regional = root/'regional'/f'islet_markers_{int(cid):04d}.tsv'
        regional_markers = pd.read_csv(regional, sep='\t')
        for marker in MARKERS:
            row = dict(case_id=cid, marker=marker, gate_log1p=np.nan, raw_threshold=np.nan,
                       comparator='>', scope='donor_all_measured_cells', validated=False,
                       evidence_status='missing_threshold_and_validation', source='', source_file='', source_sha256='',
                       regional_threshold_count=0)
            if marker in previous.index:
                old = previous.loc[marker]
                if isinstance(old, pd.DataFrame): raise ValueError('Duplicate donor/marker gates')
                supported = old.tier == 1 and not any(v in str(old.source).lower() for v in ['manual', 'provisional']) and marker != 'PDPN'
                row.update(gate_log1p=float(old.gate), raw_threshold=float(np.expm1(old.gate)),
                           validated=bool(supported), source=old.source, source_file=str(old_gates), source_sha256=legacy_hash,
                           evidence_status='existing_source_evidence' if supported else 'archived_gate_requires_review')
            if marker in hormone:
                raw = hormone[marker]
                if not np.isfinite(raw) or raw < 0: raise ValueError('Invalid hormone threshold')
                row.update(raw_threshold=raw, gate_log1p=float(np.log1p(raw)), source='qupath_explicit_threshold',
                           source_file=str(composition), source_sha256=sha256(composition))
                if not row['validated']: row['evidence_status'] = 'export_reconstructed_biological_review_pending'
            regional_values = regional_markers.loc[regional_markers.marker == marker, 'threshold']
            if len(regional_values):
                row['regional_threshold_count'] = int(regional_values.nunique())
                if not np.isfinite(row['gate_log1p']):
                    row.update(source='region_specific_export_not_a_donor_gate', source_file=str(regional),
                               source_sha256=sha256(regional), evidence_status='donor_gate_missing_region_specific_evidence_available')
            rows.append(row)
    return pd.DataFrame(rows)


def apply_reviews(registry, path):
    if path is None: return registry
    reviewed = pd.read_csv(path, dtype={'case_id': str}).fillna('')
    required = {'case_id', 'marker', 'raw_threshold', 'scope', 'decision', 'evidence_file', 'evidence_sha256', 'review_note'}
    if not required.issubset(reviewed) or reviewed.duplicated(['case_id', 'marker']).any():
        raise ValueError('Incomplete or duplicate gate review decisions')
    for r in reviewed.itertuples(index=False):
        if r.decision != 'validated': continue
        selected = (registry.case_id == r.case_id) & (registry.marker == r.marker)
        if selected.sum() != 1 or r.scope != 'donor_all_measured_cells': raise ValueError('Mixed cohort or invalid reviewed gate scope')
        evidence = Path(r.evidence_file)
        if not evidence.is_absolute(): evidence = path.parent/evidence
        if not r.review_note or not evidence.is_file() or sha256(evidence) != r.evidence_sha256:
            raise ValueError('A validated gate requires concrete matching review evidence')
        threshold = float(r.raw_threshold)
        if not np.isfinite(threshold) or threshold < 0: raise ValueError('Invalid reviewed raw gate')
        registry.loc[selected, ['raw_threshold', 'gate_log1p', 'validated', 'evidence_status', 'source', 'source_file', 'source_sha256']] = [
            threshold, float(np.log1p(threshold)), True, 'reviewed', 'explicit_gate_review', str(evidence), r.evidence_sha256]
    return registry


def gate_plot(cells, gates, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    markers = ['INS', 'GCG', 'SST', 'CD20', 'CD3e', 'CD68', 'CD163', 'CD8a', 'CD4', 'HLADR']
    fig, axes = plt.subplots(2, 5, figsize=(15, 6))
    for ax, marker in zip(axes.ravel(), markers):
        v = np.log1p(cells[marker].to_numpy(float)); v = v[np.isfinite(v)]
        ax.hist(v, bins=80, color='#657e8c', density=True)
        row = gates.loc[gates.marker == marker].iloc[0]
        if np.isfinite(row.gate_log1p):
            ax.axvline(row.gate_log1p, color='#1a7c5e' if row.validated else '#b65b08', linewidth=1.5)
        else: ax.text(.97, .93, 'Donor gate missing', ha='right', va='top', transform=ax.transAxes, fontsize=8)
        ax.set(title=marker, xlabel='log1p raw cell mean', ylabel='Density')
    fig.suptitle(f'{gates.case_id.iloc[0]}: endocrine parent-gate review\nGreen: existing source evidence; orange: archived gate awaiting review. Distributions alone do not validate gates.', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, .92))
    for suffix in ['png', 'pdf']: fig.savefig(destination.with_suffix('.'+suffix), dpi=130)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-root', required=True, type=Path)
    p.add_argument('--output-root', required=True, type=Path)
    p.add_argument('--review-root', required=True, type=Path)
    p.add_argument('--rules', required=True, type=Path)
    p.add_argument('--legacy-gates', required=True, type=Path)
    p.add_argument('--gate-review', type=Path)
    a = p.parse_args(); source = a.input_root.resolve(strict=True); out = a.output_root.resolve(strict=True)
    review = a.review_root.resolve()
    if review.exists() or (out/'trajectory/counts.csv').exists(): raise ValueError('Use fresh output artifacts; preserve previous evidence')
    review.mkdir(parents=True)
    cohort = pd.read_csv(source/'source_manifest.csv', dtype=str)
    if sha256(source/'cohort.csv') != sha256(out/'cohort.csv'): raise ValueError('Mixed source and output cohort')
    rules = read_rules(a.rules)
    registry = apply_reviews(gate_registry(source, a.legacy_gates.resolve(strict=True)), a.gate_review)
    registry.to_csv(review/'gate_registry.csv', index=False)
    template = registry.loc[~registry.validated, ['case_id', 'marker', 'raw_threshold', 'scope']].copy()
    for column in ['decision', 'evidence_file', 'evidence_sha256', 'review_note']: template[column] = ''
    template.to_csv(review/'gate_review_decisions.csv', index=False)
    counts, denominators, availability, comparisons = [], [], [], []
    features = pd.read_csv(source/'trajectory/features.csv', dtype={'case_id': str})
    for donor in cohort.itertuples(index=False):
        cid = donor.case_id; folder = review/cid; folder.mkdir()
        cells = pd.read_parquet(source/'raw'/f'{cid}-cells.parquet')
        require_cohort(cells, [cid])
        measured = cells.marker_measurements_complete.to_numpy(bool)
        logs = np.log1p(cells[MARKERS].to_numpy(dtype=float))
        gates = registry.loc[registry.case_id == cid]
        thresholds = gates.set_index('marker').gate_log1p.to_dict()
        reports, masks = audit(cells.archival_phenotype.fillna('').to_numpy(), logs, MARKERS, rules,
                               thresholds, gates.loc[gates.validated, 'marker'].tolist())
        bylabel = {r['phenotype']: r for r in reports}
        assignments = cells[['case_id', 'object_id', 'archival_phenotype', 'marker_measurements_complete']].copy()
        for label in RETAINED:
            values = pd.array(masks[label], dtype='boolean')
            values[~measured | (not bylabel[label]['validated'])] = pd.NA
            assignments[label] = values
        assignments.to_parquet(folder/'assignments.parquet', index=False)
        pd.DataFrame(reports).to_csv(folder/'retained_rule_audit.csv', index=False)
        gate_plot(cells, gates, folder/'endocrine_gate_review')
        for row in reports: availability.append(dict(case_id=cid, **row))
        cores, geoms, _ = geometry_data(donor.geometry)
        tissue = geoms['Tissue'].difference(geoms.get('Exclude', shapely.GeometryCollection()))
        union = shapely.union_all([g for _, g in cores if g.is_valid and not g.is_empty])
        x, y = cells.X_centroid.to_numpy(float)/.3774, cells.Y_centroid.to_numpy(float)/.3774
        index = shapely.STRtree(shapely.points(x, y))
        saved = pd.read_parquet(source/'analysis'/f'{cid}-memberships.parquet')
        saved_groups = {i: set(g.object_id) for i, g in saved.groupby('islet_id')}
        endo_ok = all(bylabel[label]['validated'] for label in ENDOCRINE)
        endocrine = np.logical_or.reduce([masks[label] for label in ENDOCRINE])
        missing = sorted(set().union(*(set(bylabel[t]['missing_validated_gates'].split(';')) for t in ENDOCRINE))-{''})
        for name, g in cores:
            islet = cid+'_'+name
            valid = g.is_valid and not g.is_empty and g.area > 0
            regions = {'core': g.intersection(tissue) if valid else shapely.GeometryCollection(),
                       'peri': g.buffer(20/.3774).difference(union).intersection(tissue) if valid else shapely.GeometryCollection()}
            for compartment, region in regions.items():
                members = membership(index, region, x, y); use = members[measured[members]]
                if compartment == 'core':
                    actual = set(cells.object_id.iloc[members]); expected = saved_groups.get(islet, set())
                    if actual != expected: raise ValueError(f'Core membership changed: {islet}')
                    comparisons.append(dict(islet_id=islet, case_id=cid, core_members=len(members), matches_saved_source=True))
                    reason = '' if valid and endo_ok else ('invalid_geometry' if not valid else 'gate_review_required: '+','.join(missing))
                    denominators.append(dict(islet_id=islet, case_id=cid, endocrine_count=int(endocrine[use].sum()) if valid and endo_ok else np.nan,
                        validated=bool(valid and endo_ok), reason=reason, measured_core_cells=len(use),
                        definition='unique core alpha/beta/delta cell union', rule_sha256=sha256(a.rules)))
                for label in RETAINED:
                    r = bylabel[label]; ok = valid and r['validated']
                    reason = '' if ok else ('invalid_geometry' if not valid else r['reason']+': '+r['missing_validated_gates'])
                    counts.append(dict(islet_id=islet, case_id=cid, compartment=compartment, phenotype=label,
                        count=int(masks[label][use].sum()) if ok else np.nan, validated=bool(ok), reason=reason,
                        measured_cells=len(use), valid_area_um2=region.area*.3774**2))
        print(cid, 'count artifacts built;', sum(r['validated'] for r in reports), 'of 12 retained types supported', flush=True)
    c, d, av = pd.DataFrame(counts), pd.DataFrame(denominators), pd.DataFrame(availability)
    if set(d.islet_id) != set(features.islet_id): raise ValueError('Count islets differ from frozen coordinates')
    c.to_csv(out/'trajectory/counts.csv', index=False)
    d.to_csv(out/'trajectory/endocrine_counts.csv', index=False)
    av.to_csv(out/'trajectory/phenotype_availability.csv', index=False)
    pd.DataFrame(comparisons).to_csv(review/'core_membership_verification.csv', index=False)
    write_json(out/'trajectory/marker_panel.json', dict(version=1, markers=MARKERS, measurement='per-islet median log1p core intensity'))
    write_json(out/'trajectory/count_manifest.json', dict(version=1, cohort_sha256=sha256(out/'cohort.csv'),
        rules_sha256=sha256(a.rules), gate_registry_sha256=sha256(review/'gate_registry.csv'),
        counts_sha256=sha256(out/'trajectory/counts.csv'), endocrine_counts_sha256=sha256(out/'trajectory/endocrine_counts.csv'),
        archival_labels_required=False, retained_types=RETAINED, endocrine_union=ENDOCRINE,
        wholly_unmeasured_objects_excluded=201, validated_donor_type_pairs=int(av.validated.sum()),
        supported_endocrine_donors=sorted(d.loc[d.validated, 'case_id'].unique().tolist()),
        status='available' if d.validated.all() else 'gate_review_required', source_review_root=str(review)))
    (review/'README.md').write_text('''# Retained phenotype gate review

The gate registry covers all 21 donors. Existing validated source evidence is
preserved; new exported thresholds establish computational reconstruction,
not biological validation. Regional immune thresholds are not donor-wide gates.

Each donor folder contains endocrine parent-gate plots, an archival-label
comparison, and nullable retained assignment flags. Missing assignments remain
missing. The current Ductal rule needs a reviewed positive CK19/panCK condition.

Complete gate_review_decisions.csv only for reviewed donor-wide thresholds,
with the evidence file, its SHA-256 and a review note. Rebuild into a new revision
with --gate-review pointing to that file. Parent dependencies must all pass;
new-cell assignments never require historical labels. No retired Endocrine label
or hormone-only proxy is substituted for the alpha/beta/delta union.
''')


if __name__ == '__main__': main()
