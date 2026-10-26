#!/usr/bin/env python3
"""Independent postfit diagnostics for the immune/type/size H1 extension.

No outcome models are fitted and no analysis/source files are modified. This
checks saved joins, completeness, BH arithmetic and adjusted exposure support.
Information concentration is a heuristic diagnostic, never a replacement df.
"""
from __future__ import annotations

import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.linalg as la
from statsmodels.stats.multitest import multipletests

import h1_immune_types_size_20261023 as types
import h1_immune_markers_size_20261023 as markers


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    def default(x):
        if isinstance(x, (np.integer, np.floating, np.bool_)):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, default=default, allow_nan=False) + '\n')


def bh_check(result, scope):
    rows = []
    for family, d in result.groupby('family'):
        p = d.p_cluster_t.to_numpy(float)
        expected = multipletests(np.where(np.isfinite(p), p, 1.), method='fdr_bh')[1]
        observed = d.q_bh.to_numpy(float)
        measured = np.isfinite(p)
        if not np.allclose(expected[measured], observed[measured], rtol=1e-10, atol=1e-12):
            raise ValueError(f'{scope}/{family}: BH mismatch')
        if np.isfinite(observed[~measured]).any():
            raise ValueError(f'{scope}/{family}: untested slot has a displayed q')
        rows.append(dict(scope=scope, family=family, slots=len(d), fitted=int(measured.sum()),
                         q_below_005=int(np.sum(observed[measured] < .05)), verified=True,
                         minimum_q=float(np.min(observed[measured])) if measured.any() else None))
    return rows


def info_diagnostic(d, X, detail, estimate, metadata):
    """Partial out ALL nuisance columns under saved fitted Poisson information.

    With W_i=w_i*mu_i and x the allocation column, r=x-Z(Z'WZ)^-Z'Wx.
    Donor information I_g=sum_{i in g} W_i*r_i^2. A concentration heuristic is
    (sum_g I_g)^2/sum_g I_g^2; it does not establish effective statistical df.
    """
    if list(X.columns) != detail['design_columns']:
        raise ValueError(f'Design column mismatch: {metadata}')
    a = X.to_numpy(float)
    b = np.asarray(detail['coefficients'], float)
    j = X.columns.get_loc('allocation_10pp')
    if not np.isclose(b[j], estimate['log_rr'], rtol=1e-9, atol=1e-10):
        raise ValueError('Saved coefficient does not match result row')
    y = d.n_beta_core.to_numpy(float)
    mu = np.exp(a @ b + np.log(d.core_area_um2.to_numpy(float) / 1e6))
    w = 1 / d.case_id.map(d.groupby('case_id').size()).to_numpy(float)
    w /= w.mean()
    W = w * mu
    z = np.delete(a, j, axis=1)
    x = a[:, j]
    r = x - z @ (la.pinvh((z.T * W) @ z) @ (z.T @ (W * x)))
    donors = sorted(d.case_id.unique())
    information, efficient_scores, donor_rows = [], [], []
    for donor in donors:
        mask = d.case_id.eq(donor).to_numpy()
        ig = float(np.sum(W[mask] * r[mask] ** 2))
        sg = float(np.sum(w[mask] * r[mask] * (y[mask] - mu[mask])))
        information.append(ig)
        efficient_scores.append(sg)
        donor_rows.append(dict(**metadata, case_id=donor,
                               donor_status=d.loc[mask, 'donor_status'].iloc[0],
                               objects=int(mask.sum()), core_positive_objects=int(d.loc[mask, 'type_core'].gt(0).sum()),
                               type_core_cells=int(d.loc[mask, 'type_core'].sum()),
                               beta_cells=int(d.loc[mask, 'n_beta_core'].sum()),
                               adjusted_exposure_information=ig, fitted_efficient_cluster_score=sg))
    I = np.asarray(information)
    s = np.asarray(efficient_scores)
    if I.sum() <= 0 or not np.isfinite(I).all():
        raise ValueError('No adjusted exposure information')
    shares = I / I.sum()
    s2 = s ** 2
    for row, ishare, sshare in zip(donor_rows, shares, s2 / s2.sum()):
        row['information_share'] = float(ishare)
        row['squared_fitted_score_share'] = float(sshare)
    top = int(np.argmax(I))
    score_top = int(np.argmax(s2))
    weighted_center = np.sum(W * x) / np.sum(W)
    raw_info = np.sum(W * (x - weighted_center) ** 2)
    summary = dict(**metadata, objects=len(d), donors=len(donors),
                   raw_informative_donors=int(d.groupby('case_id').type_allocation.nunique().gt(1).sum()),
                   information_concentration_donor_heuristic=float(I.sum() ** 2 / np.sum(I ** 2)),
                   maximum_information_share=float(shares[top]), maximum_information_donor=donors[top],
                   maximum_squared_fitted_score_share=float(np.max(s2) / np.sum(s2)),
                   maximum_squared_fitted_score_donor=donors[score_top],
                   squared_score_concentration_heuristic=float(np.sum(s2) ** 2 / np.sum(s2 ** 2)),
                   information_share_above_001_donors=int(np.sum(shares > .001)),
                   adjusted_to_raw_weighted_exposure_information=float(I.sum() / raw_info) if raw_info > 0 else None,
                   design_condition_number=float(np.linalg.cond(a)), rr=float(estimate['rr']),
                   saved_p=float(estimate['p_cluster_t']), saved_q=float(estimate['q_bh']),
                   covariance_or_p_values_changed=False)
    return summary, donor_rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('output/h1_immune_size_20261023'))
    p.add_argument('--features', type=Path, default=Path('output/h1_20261022/h1_object_features.csv'))
    a = p.parse_args()
    out = a.root / 'audit'
    out.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    base = types.input_data(a.features)
    table = pd.read_csv(a.root / 'object_marker_counts.csv', dtype={'islet_id': str, 'case_id': str})
    marker_results = pd.read_csv(a.root / 'marker_models/estimates.csv')
    type_results = pd.read_csv(a.root / 'types/estimates.csv')
    marker_details = json.loads((a.root / 'marker_models/fit_details.json').read_text())
    type_details = json.loads((a.root / 'types/fit_details.json').read_text())
    source_audit = json.loads((a.root / 'marker_feature_audit.json').read_text())
    dictionary = json.loads((a.root / 'data_dictionary.json').read_text())
    if not source_audit['full_cohort'] or not source_audit['all_typed_counts_match_previous_H1'] or not source_audit['all_source_cell_ids_resolved']:
        raise ValueError('Canonical source audit incomplete')
    if markers.MARKERS != dictionary['immune_associated_and_shared_state_markers'] or len(dictionary['all31_markers']) != 31:
        raise ValueError('Marker scope mismatch')
    if set(table.compartment) == {'islet', 'peri'}:
        table['compartment'] = table.compartment.replace({'islet': 'core'})
    if set(table.compartment) != {'core', 'peri'}:
        raise ValueError('Invalid compartment schema')
    if len(table) != len(base) * 2 * (31 + 12 * 2):
        raise ValueError('Marker counts table incomplete')
    if table.duplicated(['islet_id', 'compartment', 'marker', 'threshold']).any():
        raise ValueError('Duplicate object marker key')
    for population in ['all', 'immune']:
        n, measured, positive = [table[f'n_{population}{suffix}'] for suffix in ['', '_measured', '_positive']]
        if (measured < 0).any() or (measured > n).any() or (positive < 0).any() or (positive > measured).any():
            raise ValueError(f'Invalid marker counts for {population}')
        if ((n == 0) & ((measured != 0) | (positive != 0))).any():
            raise ValueError('Empty population wrongly encoded')
    bh_rows = bh_check(type_results, 'types') + bh_check(marker_results, 'markers')
    pd.DataFrame(bh_rows).to_csv(out / 'bh_verification.csv', index=False)
    completeness, summaries, donor_information = [], [], []
    for marker in markers.MARKERS:
        mf = markers.marker_frame(base, table, marker, .10)
        for size in markers.BINS:
            for cohort in ['All', 'ND', 'Aab+', 'T1D']:
                d = mf.loc[mf.equivalent_diameter.eq(size)]
                if cohort != 'All':
                    d = d.loc[d.donor_status.eq(cohort)]
                e = d.loc[d.base_eligible]
                completeness.append(dict(marker=marker, size_bin=size, cohort=cohort,
                                          source_objects=len(d), base_eligible_objects=len(e),
                                          complete_immune_marker_objects=int(e.marker_complete.sum()),
                                          incomplete_immune_marker_objects=int((~e.marker_complete).sum()),
                                          missing_core_immune_measurements=int((e.n_immune_core - e.n_immune_measured_core).sum()),
                                          missing_peri_immune_measurement_memberships=int((e.n_immune_peri - e.n_immune_measured_peri).sum()),
                                          note='Peri counts are object memberships, not unique tissue-cell counts'))
        selected = marker_results.loc[(marker_results.marker == marker) &
                                      marker_results.family.eq('conditional_broad') &
                                      marker_results.status.eq('FITTED_EXPLORATORY')]
        for _, row in selected.iterrows():
            before = mf.loc[mf.equivalent_diameter.eq(row.size_bin) & mf.base_eligible &
                            mf.marker_complete & mf.type_nearby.gt(0)].copy()
            totals = before.groupby('case_id').n_beta_core.sum()
            d = before.loc[~before.case_id.isin(totals[totals == 0].index)].copy()
            X = types.build_design(d)
            X['broad_allocation_10pp'] = d.n_immune_core / d.n_immune_nearby / .1
            key = f'conditional_broad|0.1|{marker}|{row.size_bin}|All'
            meta = dict(scope='conditional_marker', marker_or_type=marker, size_bin=row.size_bin)
            summary, donors = info_diagnostic(d, X, marker_details[key], row, meta)
            summaries.append(summary)
            donor_information.extend(donors)
    selected = type_results.loc[type_results.size_definition.eq('equivalent_diameter') &
                                type_results.cohort.eq('All') & type_results.status.eq('FITTED_EXPLORATORY')]
    for _, row in selected.iterrows():
        typed = types.type_data(base.loc[base.equivalent_diameter.eq(row.size_bin)], types.TYPES[row.immune_type])
        before = typed.loc[typed.base_eligible & typed.type_nearby.gt(0)].copy()
        totals = before.groupby('case_id').n_beta_core.sum()
        d = before.loc[~before.case_id.isin(totals[totals == 0].index)].copy()
        X = types.build_design(d)
        key = f'equivalent_diameter|{row.size_bin}|{row.immune_type}|All'
        meta = dict(scope='all_cohort_type', marker_or_type=row.immune_type, size_bin=row.size_bin)
        summary, donors = info_diagnostic(d, X, type_details[key], row, meta)
        summaries.append(summary)
        donor_information.extend(donors)
    pd.DataFrame(completeness).to_csv(out / 'marker_completeness_by_size_group.csv', index=False)
    support = pd.DataFrame(summaries)
    support.to_csv(out / 'adjusted_exposure_information_summary.csv', index=False)
    pd.DataFrame(donor_information).to_csv(out / 'adjusted_exposure_information_by_donor.csv', index=False)
    files = [a.features, a.root / 'object_marker_counts.csv', a.root / 'marker_feature_audit.json',
             a.root / 'types/estimates.csv', a.root / 'types/fit_details.json',
             a.root / 'marker_models/estimates.csv', a.root / 'marker_models/fit_details.json',
             Path(__file__), Path(types.__file__), Path(markers.__file__)]
    report = dict(recorded_utc=started, status='Independent postfit code/data/arithmetic audit; no new outcome fits or p-values',
                  objects=len(base), donors=int(base.case_id.nunique()), marker_count_rows=len(table),
                  audited_common_markers=31, model_markers=markers.MARKERS,
                  geometric_size_bins=markers.BINS, source_audit_verified=True,
                  phenotype_union_and_merge_verified=True, marker_count_bounds_verified=True,
                  bh_family_verification=bh_rows,
                  information_method='W=w*mu from saved fit; residualize allocation against every other model column using W; donor I=sum(W*r^2); heuristic=(sum I)^2/sum(I^2)',
                  information_warning='Heuristic concentration diagnostic only; not effective biological sample size, calibrated degrees of freedom, adjusted confidence intervals, a new gate, or a significance test.',
                  conditional_estimand='At fixed marker-positive/negative nearby pools and broad immune allocation, increasing marker-positive core allocation exchanges marker-positive for marker-negative cells inside and the opposite in peri. Conditional composition/segregation association, not additional infiltration or a protective/harmful effect.',
                  selection_warning='Marker/type exposed populations and zero-Beta-donor exclusions differ across fits; estimates are not direct head-to-head cell-type comparisons.',
                  size_warning='Stratified estimates do not test effect modification by size. Significance in one bin and absence/unsupportedness in another is not an interaction test.',
                  support_gate_warning='Six raw varying donors,30 core-positive objects,100 objects do not guarantee conditional exposure support or finite-cluster calibration.',
                  missingness='Complete immune-cell marker coverage is required in each object model; unknowns are not zeros. Donor descriptive means use finite measured cells; population/assay coverage must accompany comparisons.',
                  ki67='Source nuclear Ki67 is continuous within-donor immune rank, not positivity/proliferation or between-donor global change. Current curated immune selection includes Ki67; independent omission and image/control review remain necessary.',
                  inputs_and_scripts_sha256={str(p): sha(p) for p in files})
    write_json(out / 'audit_summary.json', report)
    lines = ['# Independent postfit audit — immune marker/type and size H1 extension', '',
             f'Recorded {started}. No outcome models refitted and no analysis/source inputs changed.', '',
             'Identity/union/marker completeness/count bounds and BH arithmetic passed. All36 primary marker slots remain in BH; unavailable tests contribute internal p=1 and retain missing displayed p/q. Families remain separate, as frozen.', '',
             'Adjusted exposure information was reconstructed from saved coefficients with the original donor-balanced weights. For allocation x, nuisance design Z and W=w×fitted mean, residual r=x−Z(ZᵀWZ)⁻ZᵀWx. Donor information is I_g=ΣWr². The concentration diagnostic (ΣI_g)²/ΣI_g² measures concentration only and must not replace formal degrees of freedom.', '',
             'Conditional-on-broad allocation effects are composition/segregation contrasts. They are not effects of adding immune cells, T-cell protection, LGALS3 toxicity, or proven trafficking. Model populations differ by marker/type pool and missingness. Stratified results do not establish a size interaction.', '',
             'Source nuclear Ki67 ranks remain descriptive continuous within-donor readouts. Neither ranks nor high tails prove proliferation; typing-selection dependence, localization and biological controls remain material.', '',
             '| Scope | Marker/type | Diameter bin | Donors | Raw varying donors | Information heuristic | Largest information share | Donor |',
             '|---|---|---|---:|---:|---:|---:|---|']
    for r in summaries:
        lines.append(f"| {r['scope']} | {r['marker_or_type']} | {r['size_bin']} | {r['donors']} | {r['raw_informative_donors']} | {r['information_concentration_donor_heuristic']:.2f} | {100*r['maximum_information_share']:.1f}% | {r['maximum_information_donor']} |")
    (out / 'README.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps(dict(verified=True, output=str(out), models_with_support_diagnostics=len(support),
                          smallest_information_heuristic=float(support.information_concentration_donor_heuristic.min())), indent=2))


if __name__ == '__main__':
    main()
