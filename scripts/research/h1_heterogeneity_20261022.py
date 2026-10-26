#!/usr/bin/env python3
"""Posthoc H1 group-slope and target-availability diagnostics.

Motivated by the observed positive T1D-only slope after the primary all-donor
analysis. This script does not modify, replace, or independently confirm H1.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy import stats
import statsmodels.api as sm


SCRIPT = Path(__file__).resolve()
SPEC = importlib.util.spec_from_file_location('h1_primary_20261022', SCRIPT.with_name('h1_analysis_20261022.py'))
PRIMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRIMARY)
GROUPS = ['ND', 'Aab+', 'T1D']
POSTHOC = ('Posthoc exploratory diagnostic motivated by the T1D-only reversal in the '
           'planned sensitivity analysis; primary analysis unchanged. Cross-sectional '
           'association, not target-specific killing, temporal progression, or causation.')


def write_json(path, value):
    def default(x):
        if isinstance(x, (np.integer, np.floating, np.bool_)):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value, indent=2, default=default, allow_nan=False) + '\n')


def fit_group_slopes(d):
    """Original design/weights, with one common fraction slope replaced by three."""
    zero = d.groupby('case_id').n_beta_core.sum()
    zero_donors = zero[zero.eq(0)].index.tolist()
    d = d.loc[~d.case_id.isin(zero_donors)].copy()
    if set(d.donor_status) != set(GROUPS):
        raise ValueError('All three donor groups are required for this diagnostic')
    base = PRIMARY.build_design(d, 'linear')
    allocation = base['allocation_10pp'].to_numpy(float)
    X = base.drop(columns='allocation_10pp').copy()
    names = []
    for group in GROUPS:
        col = 'allocation_10pp_' + group
        X[col] = allocation * d.donor_status.eq(group).to_numpy(float)
        names.append(col)
    a = X.to_numpy(float)
    if np.linalg.matrix_rank(a) != a.shape[1]:
        raise ValueError('Group-slope design lacks full rank')
    y = d.n_beta_core.to_numpy(float)
    offset = np.log(d.core_area_um2.to_numpy(float) / 1e6)
    w = 1 / d.case_id.map(d.groupby('case_id').size()).to_numpy(float)
    w /= w.mean()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit = sm.GLM(y, a, family=sm.families.Poisson(), offset=offset,
                     freq_weights=w).fit(maxiter=100, tol=1e-9)
    if not fit.converged or not np.isfinite(fit.params).all():
        raise ValueError('Group-slope QMLE did not converge')
    b = np.asarray(fit.params)
    mu = np.asarray(fit.fittedvalues)
    gids = d.case_id.to_numpy()
    donors = sorted(set(gids))
    G, N, P = len(donors), len(d), len(b)
    score = a * (w * (y - mu))[:, None]
    u = np.stack([score[gids == g].sum(axis=0) for g in donors])
    information = (a.T * (w * mu)) @ a
    bread = la.pinvh(information)
    correction = G / (G - 1) * (N - 1) / (N - P)
    cov = correction * bread @ (u.T @ u) @ bread
    cov = (cov + cov.T) / 2
    ix = np.array([X.columns.get_loc(name) for name in names])
    sb, scov = b[ix], cov[np.ix_(ix, ix)]
    q = stats.t.ppf(.975, G - 1)
    rows = []
    for j, group in enumerate(GROUPS):
        sub = d.loc[d.donor_status.eq(group)]
        se = float(np.sqrt(max(scov[j, j], 0)))
        rows.append(dict(donor_status=group, log_rr_10pp=sb[j], se_cluster=se,
                         rr_10pp=np.exp(sb[j]), lower_95=np.exp(sb[j] - q * se),
                         upper_95=np.exp(sb[j] + q * se), group_donors=sub.case_id.nunique(),
                         objects=len(sub), core_positive_objects=int(sub.has_core_immune.sum()),
                         ci_df=G - 1, status='posthoc exploratory, approximate cluster-t interval'))
    comparisons = [('Aab+ versus ND', np.array([-1., 1., 0.])),
                   ('T1D versus ND', np.array([-1., 0., 1.])),
                   ('T1D versus Aab+', np.array([0., -1., 1.]))]
    contrasts = []
    for name, c in comparisons:
        estimate = float(c @ sb)
        se = float(np.sqrt(max(c @ scov @ c, 0)))
        contrasts.append(dict(contrast=name, difference_log_rr_10pp=estimate,
                              se_cluster=se, ratio_of_density_ratios=np.exp(estimate),
                              lower_95=np.exp(estimate - q * se), upper_95=np.exp(estimate + q * se),
                              p_cluster_t_approx=2 * stats.t.sf(abs(estimate / se), G - 1) if se > 0 else None,
                              df=G - 1, status='posthoc exploratory, no multiplicity correction'))
    C = np.array([[-1., 1., 0.], [-1., 0., 1.]])
    cb = C @ sb
    ccov = C @ scov @ C.T
    if np.linalg.matrix_rank(ccov) != 2:
        raise ValueError('Interaction contrast covariance lacks rank two')
    wald = float(cb @ la.solve(ccov, cb, assume_a='sym'))
    F = wald / 2
    # Nested-model sanity: assigning all three slopes the primary value must
    # exactly reproduce the primary linear predictor (before optimization).
    primary_rows, primary_details, primary_fit = PRIMARY.fit_model(d, 'Primary sanity')
    common_b = primary_fit['beta']
    common_cols = primary_fit['X'].columns.tolist()
    common_params = np.array([common_b[common_cols.index(k)] if k in common_cols
                              else common_b[common_cols.index('allocation_10pp')]
                              for k in X.columns])
    common_mu = np.exp(a @ common_params + offset)
    nested_error = float(np.max(np.abs(common_mu - primary_fit['fitted'])))
    if not np.allclose(common_mu, primary_fit['fitted'], rtol=1e-10, atol=1e-10):
        raise ValueError('Nested design sanity check failed')
    if fit.llf < sm.GLM(y, base.to_numpy(float), family=sm.families.Poisson(), offset=offset,
                        freq_weights=w).loglike(common_b) - 1e-6:
        raise ValueError('Larger model has a lower fitted pseudo-log-likelihood')
    detail = dict(status=POSTHOC, objects=N, donors=donors, removed_all_zero_outcome_donors=zero_donors,
                  donor_total_weights={g: float(w[gids == g].sum()) for g in donors},
                  design_columns=X.columns.tolist(), coefficients=b.tolist(),
                  covariance_cluster=cov.tolist(), group_slope_order=GROUPS,
                  group_slope_covariance=scov.tolist(), convergence=bool(fit.converged),
                  warnings=sorted(set(str(z.message) for z in caught)),
                  nested_primary_max_mu_difference=nested_error,
                  reproduced_primary_log_rr=primary_rows[0]['log_rr'],
                  joint_equal_slopes_test=dict(wald_chi_square=wald, F_approx=F, numerator_df=2,
                                               denominator_df=G - 1, p_approx=float(stats.f.sf(F, 2, G - 1)),
                                               caveat='Approximate posthoc cluster Wald/F test; only 7,7,6 donors support the group slopes.'),
                  uncertainty_caveat='Intervals use the full model donor-cluster sandwich and t(G-1); group-specific information is limited, not 20 independent donors per group.')
    return pd.DataFrame(rows), pd.DataFrame(contrasts), detail


def census(d):
    rows, positive = [], []
    for donor, sub in d.groupby('case_id'):
        status = sub.donor_status.iloc[0]
        for beta in [False, True]:
            for intrusion in [False, True]:
                state = sub.n_beta_core.gt(0).eq(beta) & sub.has_core_immune.eq(int(intrusion))
                rows.append(dict(case_id=donor, donor_status=status, beta_present=beta,
                                 core_immune_present=intrusion, objects=int(state.sum()),
                                 donor_exposed_objects=len(sub), proportion_among_exposed=float(state.mean())))
        for intrusion in [False, True]:
            p = sub.loc[sub.n_beta_core.gt(0) & sub.has_core_immune.eq(int(intrusion))]
            positive.append(dict(case_id=donor, donor_status=status, core_immune_present=intrusion,
                                 beta_positive_objects=len(p), beta_cells=int(p.n_beta_core.sum()),
                                 beta_count_mean=float(p.n_beta_core.mean()) if len(p) else None,
                                 beta_count_median=float(p.n_beta_core.median()) if len(p) else None,
                                 beta_density_mean_mm2=float(p.beta_density_mm2.mean()) if len(p) else None,
                                 beta_density_median_mm2=float(p.beta_density_mm2.median()) if len(p) else None))
    return pd.DataFrame(rows), pd.DataFrame(positive)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--features', type=Path, default=Path('output/h1_20261022/h1_object_features.csv'))
    p.add_argument('--output', type=Path, default=Path('output/h1_20261022/heterogeneity'))
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started_utc = datetime.now(timezone.utc).isoformat()
    df = pd.read_csv(args.features, dtype={'case_id': str})
    d = PRIMARY.eligible(df)
    protocol = dict(recorded_utc=started_utc, motivation='Planned T1D-only sensitivity reversal observed after primary all-donor fitting',
                    status=POSTHOC, primary_analysis_unchanged=True,
                    recording_status='Explicitly retrospective posthoc specification; first diagnostic fit was inspected before this protocol file was persisted.',
                    population='Same PRIMARY.eligible population as the main analysis',
                    objects=len(d), donors=int(d.case_id.nunique()),
                    weights='Same inverse eligible object count per donor, normalized to mean one',
                    pooled_design='Primary donor fixed effects, class, core-area allometry, band/core ratio, nearby immune count and independent field covariates; replace common allocation slope with three group-specific allocation slopes',
                    null='Equal ND, Aab+ and T1D allocation slopes, identical weights and common nuisance covariates',
                    uncertainty='Original donor-cluster sandwich correction and t(G-1) intervals; approximate two-df cluster Wald/F; descriptive only for six T1D LODO refits',
                    followups='Per-donor four-state census Beta present/absent by core immune present/absent; selected Beta-positive object summaries are descriptive',
                    no_new_bootstrap=True)
    write_json(args.output / 'posthoc_protocol.json', protocol)
    slopes, contrasts, detail = fit_group_slopes(d)
    slopes.to_csv(args.output / 'group_specific_slopes.csv', index=False)
    contrasts.to_csv(args.output / 'group_slope_contrasts.csv', index=False)
    write_json(args.output / 'pooled_group_model.json', detail)
    states, positive = census(d)
    states.to_csv(args.output / 'four_state_census_by_donor.csv', index=False)
    positive.to_csv(args.output / 'beta_positive_object_descriptives.csv', index=False)
    t1d = d.loc[d.donor_status.eq('T1D')]
    lodo, failures = [], []
    for donor in sorted(t1d.case_id.unique()):
        try:
            rows, det, _ = PRIMARY.fit_model(t1d.loc[t1d.case_id.ne(donor)], 'Posthoc T1D LODO')
            row = rows[0]
            # Preserve approximate intervals but omit p-values: this diagnostic
            # has only five remaining donors and is used for sign/influence.
            row.pop('p_cluster_t', None)
            row.pop('t', None)
            row.update(omitted_donor=donor, status='posthoc descriptive, five remaining T1D donors')
            lodo.append(row)
        except Exception as exc:
            failures.append(dict(omitted_donor=donor, error=str(exc)))
    pd.DataFrame(lodo).to_csv(args.output / 't1d_leave_one_donor_out.csv', index=False)
    summary = dict(status=POSTHOC, primary_analysis_unchanged=True,
                   run_started_utc=started_utc, protocol_file='posthoc_protocol.json',
                   primary_population_objects=len(d), primary_population_donors=d.case_id.nunique(),
                   group_slopes=slopes.to_dict('records'),
                   joint_equal_slopes_test=detail['joint_equal_slopes_test'],
                   reproduced_primary_log_rr=detail['reproduced_primary_log_rr'],
                   nested_primary_max_mu_difference=detail['nested_primary_max_mu_difference'],
                   t1d_lodo_requested=t1d.case_id.nunique(), t1d_lodo_completed=len(lodo),
                   t1d_lodo_failed=failures,
                   t1d_lodo_rr_range=[float(min(r['rr'] for r in lodo)), float(max(r['rr'] for r in lodo))] if lodo else None,
                   four_state_census='All eligible nearby-exposed objects, including Beta-zero objects; Beta absence is an observed typed-cell state, not proof of cell death or past injury.',
                   beta_positive_descriptives='Outcome-conditioned summaries; descriptive selection, not an independent causal or inferential test.',
                   feature_sha256=hashlib.sha256(args.features.read_bytes()).hexdigest(),
                   script_sha256=hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
                   imported_primary_script_sha256=hashlib.sha256(SCRIPT.with_name('h1_analysis_20261022.py').read_bytes()).hexdigest())
    write_json(args.output / 'summary.json', summary)
    print(json.dumps(summary, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else x), flush=True)


if __name__ == '__main__':
    main()
