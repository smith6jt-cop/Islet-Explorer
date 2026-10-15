#!/usr/bin/env python3
"""Read-only verification of the frozen 15-donor trajectory and Barlow inputs.

Section 1–10 measurements describe archived H5AD DPT, not the live composition
axis (section 11). Rounded documentation values use half a final printed digit
as absolute tolerance. The known DC1/DC2 metadata discrepancy is reported
separately: preserving an old discrepancy is not evidence that it was fixed.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import importlib.metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

DONORS = {'112', '6356', '6450', '6479', '6505', '6516', '6521', '6533',
          '6534', '6538', '6548', '6549', '6550', '6551', '6563'}
STATUSES = ('ND', 'Aab+', 'T1D')
ARTIFACTS = ('data/master_results.xlsx', 'data/adata_ins_root.h5ad',
             'data/islet_explorer.h5ad')
EXPECTED = {
    'dpt_pseudotime': dict(n=5214, rho=(-.570, .143, -.206),
        donor_mean=(.3149, .3668, .3870), islet_mean=(.3040, .3634, .3883),
        donor_eta2=(.034, .090, .024)),
    'dpt_pseudotime_combined': dict(n=5214, rho=(-.511, .103, -.225),
        donor_mean=(.2317, .2649, .2853), islet_mean=(.2245, .2619, .2867),
        donor_eta2=(.024, .078, .041)),
    'dpt_pseudotime_small': dict(n=3171, rho=(-.507, .386, -.211),
        donor_mean=(.3912, .4662, .4890)),
    'dpt_pseudotime_large': dict(n=2043, rho=(-.664, .169, -.188),
        donor_mean=(.1654, .1874, .2396)),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def eta_squared(values, groups):
    values = np.asarray(values, dtype=float)
    frame = pd.DataFrame({'value': values, 'group': np.asarray(groups)})
    mean = values.mean()
    between = sum(len(g) * (g.value.mean() - mean) ** 2
                  for _, g in frame.groupby('group', observed=True))
    return float(between / np.square(values - mean).sum())


def audit(repo, baseline):
    checks = []
    evidence = dict(script_sha256=sha256(__file__), versions={name: importlib.metadata.version(name)
        for name in ('anndata', 'numpy', 'pandas', 'scipy')})

    def check(name, observed, expected, atol=None):
        passed = (bool(np.allclose(observed, expected, rtol=0, atol=atol))
                  if atol is not None else bool(observed == expected))
        checks.append(dict(check=name, passed=passed, observed=observed,
                           expected=expected, absolute_tolerance=atol))

    manifest = json.loads(baseline.read_text())
    records = {row['path']: row for row in manifest['files']}
    for name in ARTIFACTS:
        check(f'preserved_sha256:{name}', sha256(repo / name), records[name]['sha256'])
    evidence['baseline_manifest'] = str(baseline)
    evidence['baseline_manifest_sha256'] = sha256(baseline)
    evidence['barlow_aggregate_sources_sha256'] = {f'islet_analysis/islets_{region}_fixed.h5ad':
        sha256(repo / f'islet_analysis/islets_{region}_fixed.h5ad') for region in ('core', 'peri')}

    trajectory = ad.read_h5ad(repo / 'data/adata_ins_root.h5ad')
    obs = trajectory.obs.copy()
    obs['imageid'] = obs.imageid.astype(str)
    check('trajectory_shape', list(trajectory.shape), [5214, 31])
    check('original_donor_identity', sorted(obs.imageid.unique()), sorted(DONORS))
    check('unique_islet_ids', bool(obs.combined_islet_id.is_unique), True)
    params = trajectory.uns['neighbors']['params']
    for name, expected in [('use_rep', 'X_scVI_harmony_pt'), ('n_neighbors', 30),
                           ('metric', 'cosine')]:
        check('neighbor_' + name, params[name], expected)
    root = int(trajectory.uns['iroot'])
    check('root_index', root, 703)
    check('root_donor', obs.iloc[root].imageid, '6516')
    matrix = trajectory.X.toarray() if sparse.issparse(trajectory.X) else np.asarray(trajectory.X)
    marker_indices = [list(trajectory.var_names).index(m) for m in ('INS', 'GCG', 'SST')]
    check('root_INS', float(matrix[root, marker_indices[0]]), .663, .0005)
    check('first_five_diffmap_eigenvalues', trajectory.uns['diffmap_evals'][:5].tolist(),
          [1., .9399, .9268, .9142, .9060], .00005)

    for column, expected in EXPECTED.items():
        values = obs[column].to_numpy(float)
        finite = np.isfinite(values)
        frame = obs.loc[finite].copy()
        check(column + ':finite_rows', int(finite.sum()), expected['n'])
        correlations = [float(stats.spearmanr(values[finite], matrix[finite, index]).statistic)
                        for index in marker_indices]
        check(column + ':INS_GCG_SST_spearman', correlations, list(expected['rho']), .0005)
        donor_means = frame.groupby(['imageid', 'donor_status'], observed=True)[column].mean().reset_index()
        means = donor_means.groupby('donor_status', observed=True)[column].mean()
        check(column + ':donor_means_ND_Aab_T1D', [float(means[s]) for s in STATUSES],
              list(expected['donor_mean']), .00005)
        if 'islet_mean' in expected:
            means = frame.groupby('donor_status', observed=True)[column].mean()
            check(column + ':islet_means_ND_Aab_T1D', [float(means[s]) for s in STATUSES],
                  list(expected['islet_mean']), .00005)
            eta = [eta_squared(frame.loc[frame.donor_status.eq(s), column],
                               frame.loc[frame.donor_status.eq(s), 'imageid']) for s in STATUSES]
            check(column + ':within_status_donor_eta2', eta, list(expected['donor_eta2']), .0005)

    core = obs.dpt_pseudotime.to_numpy(float)
    eta_donor = eta_squared(core, obs.imageid)
    eta_status = eta_squared(core, obs.donor_status)
    check('global_donor_eta2', eta_donor, .0913, .00005)
    check('global_status_eta2', eta_status, .0445, .00005)
    check('donor_to_status_eta2_ratio', eta_donor / eta_status, 2.05, .005)
    check('core_combined_spearman', float(stats.spearmanr(core, obs.dpt_pseudotime_combined).statistic), .846, .0005)
    check('core_total_cells_spearman', float(stats.spearmanr(core, obs.total_cells).statistic), -.274, .0005)
    for threshold, expected in [(1, .111), (10, .421), (25, .155)]:
        frame = obs.loc[obs.total_cells.ge(threshold)]
        donor = frame.groupby(['imageid', 'donor_status'], observed=True).dpt_pseudotime.mean().reset_index()
        pvalue = float(stats.mannwhitneyu(donor.loc[donor.donor_status.eq('ND'), 'dpt_pseudotime'],
            donor.loc[donor.donor_status.eq('Aab+'), 'dpt_pseudotime'], alternative='less', method='exact').pvalue)
        check(f'donor_ND_Aab_exact_one_sided_MW:min_cells_{threshold}', pvalue, expected, .0005)
    check('small_stratum_mask', bool(np.array_equal(np.isfinite(obs.dpt_pseudotime_small), obs.total_cells.le(10))), True)
    check('large_stratum_mask', bool(np.array_equal(np.isfinite(obs.dpt_pseudotime_large), obs.total_cells.gt(10))), True)
    check('stratum_labels', bool(np.array_equal(obs.cell_count_stratum.astype(str),
          np.where(obs.total_cells.le(10), 'small', 'large'))), True)

    diffusion = np.asarray(trajectory.obsm['X_diffmap'])
    dc = {column: dict(matches_documented_nontrivial=bool(np.array_equal(obs[column], diffusion[:, expected])),
                      matches_preexisting_column=bool(np.array_equal(obs[column], diffusion[:, old])),
                      documented_column=expected, preexisting_column=old)
          for column, expected, old in [('DC1', 1, 0), ('DC2', 2, 1)]}
    evidence['diffusion_coordinate_metadata'] = dc
    discrepancies = []
    if not all(record['matches_documented_nontrivial'] for record in dc.values()):
        preserved = all(record['matches_preexisting_column'] for record in dc.values())
        same_hash = sha256(repo / 'data/adata_ins_root.h5ad') == records['data/adata_ins_root.h5ad']['sha256']
        discrepancies.append(dict(status='preexisting_documentation_discrepancy' if preserved and same_hash else 'unexplained_metadata_regression',
            detail='Documentation says obs DC1/DC2 use X_diffmap columns 1/2; this preserved artifact uses 0/1. The live composition trajectory does not consume these columns.',
            preserved_baseline_hash=same_hash, documented_fix_verified=False))
        check('DC_metadata_is_preserved_baseline_discrepancy', bool(preserved and same_hash), True)

    # Reconstruct the original builder's table without running its mutating main().
    aggregates = [ad.read_h5ad(repo / f'islet_analysis/islets_{region}_fixed.h5ad', backed='r')
                  for region in ('core', 'peri')]
    try:
        c, p = [a.obs[[col for col in a.obs if col.startswith('count_')]].rename(
            columns=lambda name: name.removeprefix('count_')) for a in aggregates]
        common = c.index.intersection(p.index)
        c, p = c.loc[common], p.loc[common]
        denominator = c[['Alpha cell', 'Beta cell', 'Delta cell', 'Endocrine']].sum(axis=1)
        meta = pd.DataFrame(dict(islet_id=common, case_id=aggregates[0].obs.loc[common, 'imageid'].astype(str).to_numpy(),
            donor_status=aggregates[0].obs.loc[common, 'donor_status'].astype(str).to_numpy(),
            total_cells_core=aggregates[0].obs.loc[common, 'total_cells'].astype(int).to_numpy(),
            total_cells_peri=aggregates[1].obs.loc[common, 'total_cells'].astype(int).to_numpy(),
            n_endocrine_core=denominator.to_numpy()), index=common)
        expected = pd.concat([meta, np.log1p(c.div(denominator.clip(lower=1), axis=0)).add_prefix('core_'),
            np.log1p(p.div(denominator.clip(lower=1), axis=0)).add_prefix('peri_'),
            c.add_prefix('n_core_'), p.add_prefix('n_peri_')], axis=1).reset_index(drop=True)
        actual = pd.read_csv(repo / 'data/barlow_islet_features.csv', dtype={'case_id': str}, float_precision='round_trip')
        check('barlow_shape', list(actual.shape), [5214, 90])
        check('barlow_column_identity', actual.columns.tolist(), expected.columns.tolist())
        identity = ['islet_id', 'case_id', 'donor_status']
        check('barlow_islet_order_and_metadata', bool(actual[identity].equals(expected[identity])), True)
        numeric = expected.select_dtypes(include=np.number).columns
        count_cols = [col for col in numeric if col.startswith(('n_', 'total_cells_'))]
        check('barlow_exact_raw_counts', bool(np.array_equal(actual[count_cols], expected[count_cols])), True)
        # The original builder serializes float32 log ratios using their shortest
        # decimal strings; compare that representation, not widened float32 values.
        serialized = pd.read_csv(io.StringIO(expected.to_csv(index=False)),
                                 dtype={'case_id': str}, float_precision='round_trip')
        difference = np.abs(actual[numeric].to_numpy(float) - serialized[numeric].to_numpy(float))
        check('barlow_max_absolute_numeric_error_after_builder_serialization', float(difference.max()), 0., 1e-14)
        evidence['barlow_comparison'] = dict(rows=len(actual), columns=len(actual.columns), numeric_columns=len(numeric),
            ratio_absolute_tolerance=1e-14, counts_compared_exactly=True,
            comparison='Original builder float32 calculation followed by its CSV serialization')
    finally:
        for aggregate in aggregates:
            aggregate.file.close()
    passed = all(item['passed'] for item in checks)
    return dict(passed=passed, status='restoration_checks_passed_with_preexisting_discrepancy' if passed and discrepancies
                else 'passed' if passed else 'failed', documentation_reference='000e970:app/shiny_app/docs/trajectory.md sections 1–10',
                all_documented_metadata_fixes_verified=not discrepancies, checks=checks,
                discrepancies=discrepancies, evidence=evidence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--baseline-manifest', required=True, type=Path)
    parser.add_argument('--report', type=Path, help='JSON output outside the input data directory')
    args = parser.parse_args()
    repo = args.repo_root.resolve(strict=True)
    if args.report and args.report.resolve().is_relative_to(repo / 'data'):
        parser.error('Read-only audit cannot write reports under data/')
    result = audit(repo, args.baseline_manifest.resolve(strict=True))
    text = json.dumps(result, indent=2, default=lambda value: value.item() if isinstance(value, np.generic) else str(value)) + '\n'
    if args.report:
        args.report.write_text(text)
    print(text)
    raise SystemExit(0 if result['passed'] else 1)


if __name__ == '__main__':
    main()
