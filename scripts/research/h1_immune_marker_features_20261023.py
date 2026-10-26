#!/usr/bin/env python3
"""Read-only paired-local immune marker features for the 23Oct H1 extension.

Typed immune identities and marker-positive profiles remain separate. Main
positivity is the app's exact inclusive scaled>=0.10 cut;0.05/0.20 are threshold
sensitivities for12 immune-associated/shared-state markers. Donor protein means
deduplicate persistent cells; object profiles deliberately retain peri sharing.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import io
import json
from pathlib import Path
import platform
import sys
import traceback

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import h1_features_20261022 as hf

MARKERS = ['CD31', 'CD8a', 'CD3e', 'SMA', 'CD4', 'CD34', 'HLADR', 'PDPN', 'panCK', 'ECAD', 'CD163',
           'SST', 'ColIV', 'VIM', 'CD20', 'LGALS3', 'B3TUBB', 'GCG', 'KRT14', 'GAP43', 'CD35', 'CHGA',
           'PGP9.5', 'INS', 'CD44', 'CD45', 'NaKATPase', 'BCatenin', 'CD68', 'BActin', 'CK19']
IMMUNE_MARKERS = ['CD45', 'CD3e', 'CD8a', 'CD4', 'CD20', 'CD68', 'CD163', 'HLADR', 'CD35', 'LGALS3', 'CD44', 'VIM']
IMMUNE_TYPES = ['CD8a Tcell', 'CD4 Tcell', 'T cell', 'B cell', 'Macrophage', 'Immune']
POPULATIONS = ['ALL', 'IMMUNE_UNION', *IMMUNE_TYPES]
SIZE_BINS = ['<50', '50-<100', '>=100']
CUTS = [(0.10, np.arange(len(MARKERS))), (0.05, np.array([MARKERS.index(m) for m in IMMUNE_MARKERS])),
        (0.20, np.array([MARKERS.index(m) for m in IMMUNE_MARKERS]))]


def size_bin(diameter):
    return '<50' if diameter < 50 else ('50-<100' if diameter < 100 else '>=100')


def population_statistics(d):
    """Finite-cell sufficient statistics; empty pools differ from no assay."""
    X = d[MARKERS].to_numpy(float)
    labels = d.phenotype.to_numpy()
    result = {}
    for population in POPULATIONS:
        mask = np.ones(len(d), bool) if population == 'ALL' else (
            np.isin(labels, IMMUNE_TYPES) if population == 'IMMUNE_UNION' else labels == population)
        values = X[mask]
        measured = np.isfinite(values)
        n = measured.sum(axis=0).astype(np.int64)
        clean = np.where(measured, values, 0)
        sums = clean.sum(axis=0)
        sums_sq = np.square(clean).sum(axis=0)
        sums[n == 0], sums_sq[n == 0] = np.nan, np.nan
        means = np.divide(sums, n, out=np.full(len(MARKERS), np.nan), where=n > 0)
        positives = {}
        for cut, _ in CUTS:
            p = (measured & (values >= cut)).sum(axis=0).astype(float)
            if len(values):
                p[n == 0] = np.nan
            positives[cut] = p
        result[population] = dict(n_cells=int(mask.sum()), n_measured=n, sum=sums,
                                  sum_sq=sums_sq, mean=means, positives=positives)
    return result


def append_profiles(rows, base, stats, cut, indices, retain_empty=False):
    for population, s in stats.items():
        if s['n_cells'] == 0 and not retain_empty:
            continue
        for j in indices:
            measured = int(s['n_measured'][j])
            positive = float(s['positives'][cut][j])
            rows.append(dict(base, phenotype=population, marker=MARKERS[j], threshold=cut,
                             n_cells=s['n_cells'], n_measured=measured, n_missing=s['n_cells'] - measured,
                             n_positive=positive, positive_fraction=positive / measured if measured else np.nan,
                             sum=float(s['sum'][j]), sum_sq=float(s['sum_sq'][j]), mean=float(s['mean'][j])))


def donor_profiles(memberships, case, status, diameter=None):
    d = memberships if diameter is None else memberships.loc[memberships.diameter_bin.eq(diameter)]
    d = d.drop_duplicates(['object_id', 'cell_region'])
    rows = []
    for compartment in ['core', 'peri']:
        stats = population_statistics(d.loc[d.cell_region.eq(compartment)])
        base = dict(case_id=case, donor_status=status, compartment=compartment,
                    diameter_bin='ALL' if diameter is None else diameter,
                    aggregation_unit='unique_cell_within_donor_compartment' if diameter is None else
                    'unique_cell_within_donor_compartment_and_size_bin')
        append_profiles(rows, base, stats, .10, np.arange(len(MARKERS)), retain_empty=True)
    return rows


def donor_extract(root_str, out_str, case, objects, prior_hashes):
    root, out = Path(root_str), Path(out_str)
    source_path = root / 'cell_types' / f'{case}.parquet'
    source = pd.read_parquet(source_path, columns=['object_id', 'case_id', 'phenotype', 'cell_region', 'islet_name', 'x_um', 'y_um'])
    if source.object_id.duplicated().any() or set(source.case_id.astype(str)) != {str(case)}:
        raise ValueError(f'{case}: duplicate persistent identity or invalid donor')
    index = pd.MultiIndex.from_frame(source[['x_um', 'y_um', 'phenotype']])
    if not index.is_unique:
        raise ValueError(f'{case}: coordinates+phenotype cannot uniquely identify source cells')
    lookup = pd.Series(source.object_id.to_numpy(), index=index)
    source_by_id = source.set_index('object_id', drop=False)
    primary, sensitivity, compact, denominators, donor_means, memberships = [], [], [], [], [], []
    hashes = {str(source_path): hf.sha256(source_path)}
    marker_missing = np.zeros(len(MARKERS), dtype=np.int64)
    n_memberships = 0
    for rec in objects.to_dict('records'):
        path = root / 'cells' / f"{rec['islet_id']}.csv"
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if prior_hashes.get(str(path)) != digest:
            raise ValueError(f'{rec["islet_id"]}: drilldown changed since H1 audit')
        hashes[str(path)] = digest
        d = pd.read_csv(io.BytesIO(raw), float_precision='round_trip',
                        usecols=['X_centroid', 'Y_centroid', 'phenotype', 'cell_region', *MARKERS])
        if not set(d.cell_region).issubset({'core', 'peri'}) or not set(d.phenotype).issubset(hf.LABELS):
            raise ValueError(f'{path}: invalid region or phenotype')
        values = d[MARKERS].to_numpy(float)
        finite = np.isfinite(values)
        if ((values[finite] < 0) | (values[finite] > 1)).any():
            raise ValueError(f'{path}: marker is outside exact app[0,1] scale')
        marker_missing += (~finite).sum(axis=0)
        n_memberships += len(d)
        # Explicit coordinate+phenotype key lookup, preserving persistent IDs;
        # no CSV/Parquet row-order correspondence is assumed.
        key = pd.MultiIndex.from_frame(d[['X_centroid', 'Y_centroid', 'phenotype']].rename(
            columns={'X_centroid': 'x_um', 'Y_centroid': 'y_um'}))
        ids = lookup.reindex(key)
        if ids.isna().any() or ids.duplicated().any():
            raise ValueError(f'{path}: source cell identity lookup missing/ambiguous')
        d['object_id'] = ids.to_numpy()
        core = d.cell_region.eq('core').to_numpy()
        core_source = source_by_id.loc[d.loc[core, 'object_id']]
        if not core_source.cell_region.eq('core').all() or not core_source.islet_name.eq(rec['base_islet_id']).all():
            raise ValueError(f'{path}: persistent core owner mismatch')
        if source_by_id.loc[d.loc[~core, 'object_id'], 'cell_region'].eq('core').any():
            raise ValueError(f'{path}: peri table includes another endocrine core')
        diameter = float(2 * np.sqrt(rec['core_area_um2'] / np.pi))
        bin_name = size_bin(diameter)
        d['diameter_bin'] = bin_name
        memberships.append(d[['object_id', 'phenotype', 'cell_region', 'diameter_bin', *MARKERS]])
        base = dict(islet_id=rec['islet_id'], case_id=case, donor_status=rec['donor_status'],
                    object_class=rec['object_class'], core_area_um2=float(rec['core_area_um2']),
                    equivalent_diameter_um=diameter, diameter_bin=bin_name,
                    n_endocrine_core=int(rec['n_endocrine_core']), n_cells_core=int(rec['n_cells_core']))
        for compartment in ['core', 'peri']:
            region = d.loc[d.cell_region.eq(compartment)]
            actual = region.phenotype.value_counts().reindex(hf.LABELS, fill_value=0)
            for label in hf.LABELS:
                if actual[label] != rec[f'n_{compartment}_{hf.token(label)}']:
                    raise ValueError(f'{path}: {compartment}/{label} fails prior H1 count audit')
            stats = population_statistics(region)
            subbase = dict(base, compartment=compartment)
            for population, s in stats.items():
                denominators.append(dict(subbase, phenotype=population, n_cells=s['n_cells'],
                                         population_present=s['n_cells'] > 0))
            for cut, ix in CUTS:
                append_profiles(primary if cut == .10 else sensitivity, subbase, stats, cut, ix)
                all_s, immune_s = stats['ALL'], stats['IMMUNE_UNION']
                for j in ix:
                    compact.append(dict(subbase, marker=MARKERS[j], threshold=cut,
                                        n_all=all_s['n_cells'], n_all_measured=int(all_s['n_measured'][j]),
                                        n_all_positive=float(all_s['positives'][cut][j]),
                                        n_immune=immune_s['n_cells'], n_immune_measured=int(immune_s['n_measured'][j]),
                                        n_immune_positive=float(immune_s['positives'][cut][j])))
    memberships = pd.concat(memberships, ignore_index=True)
    duplicated = memberships.loc[memberships.duplicated(['object_id', 'cell_region'], keep=False)]
    if len(duplicated):
        unique_values = duplicated.groupby(['object_id', 'cell_region'], sort=False)[MARKERS].nunique(dropna=False)
        if unique_values.gt(1).any().any():
            raise ValueError(f'{case}: repeated peri-cell protein values disagree')
    unique_local = memberships.drop_duplicates(['object_id', 'cell_region'])
    if unique_local.object_id.duplicated().any():
        raise ValueError(f'{case}: same cell enters both core and peri')
    donor_means.extend(donor_profiles(memberships, case, str(objects.donor_status.iloc[0])))
    size_means = []
    for diameter in SIZE_BINS:
        size_means.extend(donor_profiles(memberships, case, str(objects.donor_status.iloc[0]), diameter))
    folder = out / 'marker_features_by_donor' / str(case)
    folder.mkdir(parents=True, exist_ok=True)
    tables = {'marker_profiles': primary, 'marker_profiles_sensitivity': sensitivity,
              'object_marker_counts': compact, 'population_denominators': denominators,
              'donor_type_compartment_markers': donor_means, 'donor_size_type_compartment_markers': size_means}
    for name, rows in tables.items():
        pd.DataFrame(rows).to_parquet(folder / f'{name}.parquet', index=False)
    audit = dict(case_id=case, donor_status=str(objects.donor_status.iloc[0]), objects=len(objects),
                 paired_cell_memberships=n_memberships, paired_unique_cells=len(unique_local),
                 paired_unique_core_cells=int(unique_local.cell_region.eq('core').sum()),
                 paired_unique_peri_cells=int(unique_local.cell_region.eq('peri').sum()),
                 typed_counts_match_previous_H1=True, source_cell_ids_resolved=True,
                 marker_missing_memberships=dict(zip(MARKERS, marker_missing.tolist())),
                 rows={k: len(v) for k, v in tables.items()}, source_sha256=hashes,
                 output_sha256={str(folder / f'{name}.parquet'): hf.sha256(folder / f'{name}.parquet') for name in tables})
    hf.write_json(folder / 'audit.json', audit)
    print(f'{case}: {len(objects)} objects/{n_memberships} memberships; persistent IDs and all typed counts verified', flush=True)
    return audit


def combine_parquet(paths, target):
    writer = None
    n = 0
    try:
        for p in paths:
            table = pq.read_table(p)
            if writer is None:
                writer = pq.ParquetWriter(target, table.schema, compression='zstd', use_dictionary=True)
            writer.write_table(table)
            n += len(table)
    finally:
        if writer is not None:
            writer.close()
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--h1-features', type=Path, required=True)
    parser.add_argument('--h1-audit', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--cases', nargs='+', type=int, help='Pilot subset; final canonical run must omit')
    parser.add_argument('--checkpoint-only', action='store_true')
    args = parser.parse_args()
    root, out = args.data_root.resolve(strict=True), args.output_root.resolve()
    if out == root or out.is_relative_to(root):
        raise ValueError('Read-only source data: output must be outside served root')
    out.mkdir(parents=True, exist_ok=True)
    objects = pd.read_csv(args.h1_features)
    old_audit = json.loads(args.h1_audit.read_text())
    feature_hash = hf.sha256(args.h1_features)
    recorded_feature_hash = old_audit['output_sha256'].get(str(args.h1_features.resolve()))
    if old_audit['data_root'] != str(root) or feature_hash != recorded_feature_hash:
        raise ValueError('Old H1 features/source provenance mismatch')
    donors = objects[['case_id', 'donor_status']].drop_duplicates()
    if len(donors) != 20 or 78 in set(donors.case_id) or objects.islet_id.duplicated().any():
        raise ValueError('Expected unchanged20-donor H1 cohort, excluding78')
    if donors.donor_status.value_counts().to_dict() != {'ND': 7, 'Aab+': 7, 'T1D': 6}:
        raise ValueError('Unexpected donor groups')
    cases = sorted(set(args.cases) if args.cases else set(donors.case_id))
    if not set(cases).issubset(set(donors.case_id)):
        raise ValueError('Absent/excluded requested case')
    audits, failures = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(donor_extract, str(root), str(out), case, objects.loc[objects.case_id.eq(case)],
                               old_audit['source_sha256']): case for case in cases}
        for future in as_completed(futures):
            case = futures[future]
            try:
                audits.append(future.result())
            except Exception as exc:
                failures.append(dict(case_id=case, error=str(exc), traceback=traceback.format_exc()))
                hf.write_json(out / 'FAILED_marker_extraction.json', failures)
                print(f'{case}: FAILED {exc}', file=sys.stderr, flush=True)
    if failures:
        raise RuntimeError('Marker extraction failed; inspect FAILED_marker_extraction.json')
    if args.checkpoint_only:
        return
    folder = out / 'marker_features_by_donor'
    names = ['marker_profiles', 'marker_profiles_sensitivity', 'object_marker_counts', 'population_denominators',
             'donor_type_compartment_markers', 'donor_size_type_compartment_markers']
    row_counts = {}
    for name in names:
        row_counts[name] = combine_parquet([folder / str(case) / f'{name}.parquet' for case in cases], out / f'{name}.parquet')
        if name not in {'marker_profiles', 'marker_profiles_sensitivity'}:
            pd.read_parquet(out / f'{name}.parquet').to_csv(out / f'{name}.csv', index=False)
    denominators = pd.read_parquet(out / 'population_denominators.parquet')
    expected_objects = objects.loc[objects.case_id.isin(cases)]
    if len(denominators) != len(expected_objects) * 2 * len(POPULATIONS):
        raise ValueError('Incomplete population denominators')
    profiles = pd.read_parquet(out / 'marker_profiles.parquet')
    coverage = profiles.loc[profiles.phenotype.eq('ALL')].groupby(['marker', 'compartment']).agg(
        objects=('islet_id', 'size'), n_cell_memberships=('n_cells', 'sum'), n_measured=('n_measured', 'sum'),
        n_missing=('n_missing', 'sum'), n_positive=('n_positive', lambda x: x.sum(min_count=1))).reset_index()
    coverage.to_csv(out / 'marker_coverage.csv', index=False)
    hashes = {str(args.h1_features.resolve()): feature_hash, str(args.h1_audit.resolve()): hf.sha256(args.h1_audit)}
    for audit in audits:
        hashes.update(audit['source_sha256'])
    dictionary = dict(keys=['islet_id', 'case_id', 'donor_status', 'object_class', 'compartment', 'phenotype', 'marker', 'threshold'],
                      populations=POPULATIONS, curated_immune_union=IMMUNE_TYPES,
                      all31_markers=MARKERS, immune_associated_and_shared_state_markers=IMMUNE_MARKERS,
                      scale='compartment REDSEA/UniFORM-aligned/cohort floor-top scaled and clipped[0,1]; scaled cell mean, not protein mass',
                      positivity='inclusive >=threshold, exactly0.10 primary for all31;0.05/0.20 sensitivity for12markers',
                      missingness='nonempty populations with zero finite measurements have n_positive,sum,sum_sq,mean missing; empty populations have n_cells0 in complete denominator table and are omitted from object marker profiles',
                      denominator='n_measured is positivity/mean denominator; n_cells includes unmeasured cells; n_positive0 for a truly empty immune pool is an empty-cell count, not a negative assay',
                      sum='sum of finite scaled per-cell means, not integrated intensity; sum_sq is sum of their squares',
                      size='diameter2*sqrt(core_area_um2/pi), bins<50,50-<100,>=100µm; n_endocrine_core/n_cells_core retain alternative size information',
                      object_aggregation='core/peri memberships preserve existing per-object shared peri bands',
                      donor_aggregation='persistent object_id deduplication within donor+compartment; size companion deduplicates within size bin, and the same peri cell can occur in different bins',
                      identity_linkage='explicit source x_um+y_um+phenotype key resolves unique persistent object_id; no row-order joins',
                      biological_limits='marker-positive populations overlap and do not define new immune identities; CD44/VIM/LGALS3/HLADR can be shared with nonimmune cells; no activation/function or Ki67 proliferation call inferred')
    hf.write_json(out / 'data_dictionary.json', dictionary)
    hf.write_json(out / 'marker_feature_audit.json', dict(full_cohort=not bool(args.cases), data_root=str(root),
                  cases=cases, objects=len(expected_objects), rows=row_counts,
                  paired_cell_memberships=sum(a['paired_cell_memberships'] for a in audits),
                  paired_unique_cells=sum(a['paired_unique_cells'] for a in audits),
                  all_typed_counts_match_previous_H1=all(a['typed_counts_match_previous_H1'] for a in audits),
                  all_source_cell_ids_resolved=all(a['source_cell_ids_resolved'] for a in audits),
                  source_sha256=hashes, script_sha256=hf.sha256(Path(__file__)),
                  helper_script_sha256=hf.sha256(Path(hf.__file__)), argv=sys.argv,
                  software={'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__, 'pyarrow': pa.__version__},
                  output_sha256={str(out / f'{name}.parquet'): hf.sha256(out / f'{name}.parquet') for name in names}))
    print(f'Canonical marker tables complete: {len(expected_objects)} objects/{len(cases)} donors; {out}', flush=True)


if __name__ == '__main__':
    main()
