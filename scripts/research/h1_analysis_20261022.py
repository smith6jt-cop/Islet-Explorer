#!/usr/bin/env python3
"""Read-only H1 analysis of the accepted uniform 20-donor data.

Run with the project's scvi-env Python. Inputs are independently verified by
h1_features_20261022.py. Outputs never enter the app's data root.
"""
from __future__ import annotations

import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/islet-h1-matplotlib')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import argparse
import hashlib
import json
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy import stats
import statsmodels.api as sm

SEED = 20261022
STATUS = ['ND', 'Aab+', 'T1D']
COLORS = {'ND': '#4477AA', 'Aab+': '#CC6633', 'T1D': '#228833'}
ADDITIONS = {'6090', '6232', '6520', '6539', '6558'}
OUT = Path('output/h1_20261022')


def write_json(path, obj):
    def default(x):
        if isinstance(x, (np.integer, np.floating, np.bool_)):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(obj, indent=2, default=default) + '\n')


def eligible(df, field='20_100', min_field_area=1000, exclusive=False):
    d = df.copy()
    d['case_id'] = d.case_id.astype(str)
    if '78' in set(d.case_id):
        raise ValueError('Excluded donor78 in input')
    if exclusive:
        d['n_immune_nearby'] = d.n_immune_nearby_exclusive
        d['immune_localization_fraction'] = d.immune_localization_fraction_exclusive
    d['field_density'] = d[f'immune_density_field{field}_mm2']
    good = ((d.n_immune_nearby > 0) & (d.core_area_um2 > 0)
            & (d.peri_area_um2 > 0)
            & (d[f'field{field}_area_um2'] >= min_field_area)
            & np.isfinite(d.field_density)
            & np.isfinite(d.immune_localization_fraction))
    edge_columns = {'core_complete', 'neighborhood_complete_20', 'neighborhood_complete_100',
                    'neighborhood_complete_250'}
    if edge_columns.intersection(d.columns):
        if not edge_columns.issubset(d.columns) or d[list(edge_columns)].isna().any().any():
            raise ValueError('Incomplete acquisition eligibility metadata')
        for column in edge_columns:
            if not d[column].isin([True, False]).all():
                raise ValueError(f'Invalid boolean eligibility: {column}')
        radius = int(field.split('_')[1])
        good &= d.core_complete & d.neighborhood_complete_20 & d[f'neighborhood_complete_{radius}']
    d = d.loc[good].copy()
    d['beta_density_mm2'] = d.n_beta_core / d.core_area_um2 * 1e6
    d['has_core_immune'] = (d.n_immune_core > 0).astype(int)
    d['nonbeta_endocrine'] = d.n_endocrine_core - d.n_beta_core
    return d


def build_design(d, exposure='linear'):
    # Fractions are in 10-percentage-point units, so beta is directly log(RR10pp).
    X = pd.DataFrame(index=d.index)
    X['intercept'] = 1.0
    for donor in sorted(d.case_id.unique())[1:]:
        X['donor_' + donor] = (d.case_id == donor).astype(float)
    X['class_cluster'] = (d.object_class != 'Islet').astype(float)
    nuisance = {
        'log_core_area': np.log(d.core_area_um2),
        'log_band_core_area_ratio': np.log(d.peri_area_um2 / d.core_area_um2),
        'log1p_nearby_immune': np.log1p(d.n_immune_nearby),
        'log1p_independent_field': np.log1p(d.field_density),
    }
    for k, v in nuisance.items():
        sd = v.std(ddof=0)
        if sd > 1e-12:
            X[k] = (v - v.mean()) / sd
    if exposure == 'linear':
        X['allocation_10pp'] = d.immune_localization_fraction / .10
    elif exposure == 'two_part':
        median = d.loc[d.has_core_immune == 1, 'immune_localization_fraction'].median()
        X['any_core_immune'] = d.has_core_immune
        X['positive_allocation_10pp'] = d.has_core_immune * (
            d.immune_localization_fraction - median) / .10
    else:
        raise ValueError(exposure)
    # Class is constant in Islets-only sensitivity. Drop exact constant columns
    # except the intercept, then drop numerical dependencies transparently.
    X = X.loc[:, [(k == 'intercept') or (X[k].std(ddof=0) > 1e-12) for k in X]]
    rank = np.linalg.matrix_rank(X.to_numpy(float))
    if rank != X.shape[1]:
        raise ValueError(f'Design not identifiable: rank {rank}/{X.shape[1]}')
    return X.astype(float)


def fit_model(data, name, outcome='n_beta_core', exposure='linear', balanced=True,
              inference=True, wild=False):
    d = data.copy()
    zero = d.groupby('case_id')[outcome].sum()
    zero_donors = zero[zero == 0].index.tolist()
    d = d.loc[~d.case_id.isin(zero_donors)].copy()
    if len(d) == 0 or d.case_id.nunique() < 2:
        raise ValueError('Insufficient nonzero-outcome donors')
    X = build_design(d, exposure)
    if exposure == 'linear' and 'allocation_10pp' not in X:
        raise ValueError('Allocation has no variation')
    y = d[outcome].to_numpy(float)
    offset = np.log(d.core_area_um2.to_numpy(float) / 1e6)
    if balanced:
        sizes = d.case_id.map(d.groupby('case_id').size()).to_numpy(float)
        w = 1 / sizes
        w = w / w.mean()
    else:
        w = np.ones(len(d))
    a = X.to_numpy(float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit = sm.GLM(y, a, family=sm.families.Poisson(), offset=offset,
                     freq_weights=w).fit(maxiter=100, tol=1e-9)
    if not fit.converged or not np.all(np.isfinite(fit.params)):
        raise ValueError('Poisson mean fit did not converge')
    b = np.asarray(fit.params)
    mu = np.asarray(fit.fittedvalues)
    groups = d.case_id.to_numpy()
    donors = sorted(set(groups))
    G, N, P = len(donors), len(d), a.shape[1]
    score = a * (w * (y - mu))[:, None]
    u = np.stack([score[groups == g].sum(axis=0) for g in donors])
    bread = la.pinvh((a.T * (w * mu)) @ a)
    correction = G / (G - 1) * (N - 1) / (N - P)
    cov = correction * bread @ (u.T @ u) @ bread
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    q = stats.t.ppf(.975, G - 1)
    rows = []
    for k in ['allocation_10pp', 'any_core_immune', 'positive_allocation_10pp']:
        if k not in X:
            continue
        j = X.columns.get_loc(k)
        t = b[j] / se[j] if se[j] > 0 else np.nan
        rows.append(dict(model=name, outcome=outcome, term=k, log_rr=b[j], se_cluster=se[j],
                         rr=np.exp(b[j]), lower=np.exp(b[j]-q*se[j]), upper=np.exp(b[j]+q*se[j]),
                         p_cluster_t=2*stats.t.sf(abs(t), G-1), t=t, donors=G,
                         objects=N, core_positive=int(d.has_core_immune.sum()),
                         informative_donors=int(d.groupby('case_id').immune_localization_fraction.nunique().gt(1).sum()),
                         donor_weights='equal object-weight per donor' if balanced else 'unweighted objects'))
    details = dict(name=name, outcome=outcome, objects=N, donors=donors,
                   all_zero_outcome_donors=zero_donors,
                   donor_total_weights={g: float(w[groups == g].sum()) for g in donors},
                   design_columns=X.columns.tolist(), coefficients=b.tolist(),
                   covariance_cluster=cov.tolist(), convergence=bool(fit.converged),
                   warnings=sorted(set(str(z.message) for z in caught)),
                   pearson_dispersion=float(np.sum((y-mu)**2 / np.maximum(mu, 1e-12)) / (N-P)),
                   allocation_positive_median=float(d.loc[d.has_core_immune == 1, 'immune_localization_fraction'].median()))
    if wild and exposure == 'linear':
        j = X.columns.get_loc('allocation_10pp')
        z = a[:, j]
        nuisance = np.delete(a, j, axis=1)
        restricted = sm.GLM(y, nuisance, family=sm.families.Poisson(), offset=offset,
                            freq_weights=w).fit(maxiter=100, tol=1e-9)
        mu0 = np.asarray(restricted.fittedvalues)
        if not restricted.converged or not np.all(np.isfinite(mu0)):
            raise ValueError('Restricted score fit did not converge')
        info0 = (nuisance.T * (w * mu0)) @ nuisance
        gamma = la.pinvh(info0) @ (nuisance.T @ (w*mu0*z))
        zres = z - nuisance @ gamma
        scores = np.array([np.sum(w[groups == g] * zres[groups == g] *
                                  (y[groups == g] - mu0[groups == g])) for g in donors])
        den = np.sqrt(np.sum(scores**2))
        if den <= 0 or not np.isfinite(den):
            raise ValueError('Wild score variance is unavailable')
        observed = scores.sum() / den
        rng = np.random.default_rng(SEED)
        simulated = rng.choice([-1., 1.], size=(9999, G)) @ scores / den
        p = (1 + np.sum(np.abs(simulated) >= abs(observed))) / 10000
        details['wild_restricted_score'] = dict(draws=9999, statistic=observed, p_two_sided=p,
            seed=SEED, cluster_scores=dict(zip(donors, scores.tolist())),
            interpretation='Approximate donor score-symmetry sensitivity, not exact randomization')
    return rows, details, dict(data=d, X=X, beta=b, covariance=cov, weights=w, fitted=mu)


def boot_one(args):
    d, i = args
    rng = np.random.default_rng(SEED + 1000 + i)
    donors = sorted(d.case_id.unique())
    draw = rng.choice(donors, size=len(donors), replace=True)
    pieces = []
    for j, g in enumerate(draw):
        sub = d.loc[d.case_id == g].copy()
        sub['case_id'] = f'bootstrap_{j:02d}'
        pieces.append(sub)
    try:
        rows, details, _ = fit_model(pd.concat(pieces, ignore_index=True), 'bootstrap', inference=False)
        return dict(draw=i, log_rr=rows[0]['log_rr'], error=None)
    except Exception as exc:
        return dict(draw=i, log_rr=np.nan, error=str(exc))


def support_table(df, d):
    rows = []
    for g, allsub in df.groupby('case_id'):
        s = d.loc[d.case_id == str(g)]
        rows.append(dict(case_id=str(g), donor_status=allsub.donor_status.iloc[0],
            paired_objects=len(allsub), exposed_objects=int((allsub.n_immune_nearby > 0).sum()),
            primary_objects=len(s), core_positive=int((s.n_immune_core > 0).sum()),
            core_immune_cells=int(s.n_immune_core.sum()), nearby_immune_memberships=int(s.n_immune_nearby.sum()),
            no_exposure_objects=int((allsub.n_immune_nearby == 0).sum()),
            allocation_unique=int(s.immune_localization_fraction.nunique()),
            allocation_max=float(s.immune_localization_fraction.max()) if len(s) else np.nan,
            beta_cells=int(s.n_beta_core.sum()), beta_zero_objects=int((s.n_beta_core == 0).sum())))
    s = pd.DataFrame(rows)
    s['status_order'] = s.donor_status.map(dict(zip(STATUS, range(3))))
    return s.sort_values(['status_order', 'case_id']).drop(columns='status_order')


def save_figure(fig, name, out):
    fig.savefig(out / (name + '.png'), dpi=200, facecolor='white', bbox_inches='tight')
    fig.savefig(out / (name + '.svg'), facecolor='white', bbox_inches='tight')
    plt.close(fig)


def figures(df, d, support, result, details, fit, lodo, out):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titlesize': 12, 'axes.labelsize': 11})
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), gridspec_kw={'width_ratios': [1.35, 1]})
    y = np.arange(len(support))
    labels = [f'{r.case_id}  {r.donor_status}' for r in support.itertuples()]
    axes[0].barh(y, support.primary_objects, color='#D9DEE5', label='No core immune cells')
    axes[0].barh(y, support.core_positive, color=[COLORS[x] for x in support.donor_status], label='Core immune cells present')
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel('Eligible immune-exposed objects')
    axes[0].set_title('A  Where the information comes from', loc='left', fontweight='bold')
    axes[0].legend(loc='lower right', fontsize=9)
    for j, r in enumerate(support.itertuples()):
        axes[0].text(r.primary_objects + 4, j, f'{r.core_positive}/{r.primary_objects}', va='center', fontsize=8)
    axes[0].set_xlim(0, support.primary_objects.max() * 1.26)
    axes[1].barh(y, support.core_immune_cells, color=[COLORS[x] for x in support.donor_status])
    axes[1].set_yticks(y, labels)
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Broad immune cells inside eligible cores')
    axes[1].set_title('B  Core immune events by donor', loc='left', fontweight='bold')
    for j, r in enumerate(support.itertuples()):
        axes[1].text(r.core_immune_cells + 1, j, str(r.core_immune_cells), va='center', fontsize=8)
    axes[1].set_xlim(0, max(support.core_immune_cells.max() * 1.18, 5))
    fig.suptitle('H1 data support: donors, objects, and sparse core immune cells', fontsize=15, fontweight='bold')
    fig.text(.5, .015, 'All counts use the analysis eligibility rules. Nearby immune cells absent: analysed separately. '
             'Bars count objects; donors remain the independent units.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .04, 1, .95])
    save_figure(fig, 'figure1_donor_support', out)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), gridspec_kw={'width_ratios': [1.05, 1]})
    # Observed donor summaries, not a pooled object-level significance test.
    points = []
    for g, sub in d.groupby('case_id'):
        a, b = sub.loc[sub.has_core_immune == 0], sub.loc[sub.has_core_immune == 1]
        if len(a) and len(b):
            status = sub.donor_status.iloc[0]
            vals = [a.beta_density_mm2.mean() / 1000, b.beta_density_mm2.mean() / 1000]
            axes[0].plot([0, 1], vals, color=COLORS[status], alpha=.55, lw=1.2)
            axes[0].scatter([0, 1], vals, color=COLORS[status], s=25, zorder=3)
            points.append(dict(case_id=g, donor_status=status, mean_density_no_core=vals[0], mean_density_core=vals[1]))
    for status in STATUS:
        axes[0].plot([], [], color=COLORS[status], marker='o', label=status)
    axes[0].set_xticks([0, 1], ['Immune cells around only', 'Some immune cells inside'])
    axes[0].set_ylabel('Mean Beta density (cells per 1,000 µm²)')
    axes[0].set_title('A  Observed within-donor comparison', loc='left', fontweight='bold')
    axes[0].legend(title='Donor group', fontsize=9)
    pd.DataFrame(points).to_csv(out / 'observed_donor_intrusion_comparison.csv', index=False)
    # The relative curve holds all model covariates and donor intercepts fixed.
    r = result.loc[result.model == 'Primary'].iloc[0]
    xmax = min(.5, d.loc[d.has_core_immune == 1, 'immune_localization_fraction'].quantile(.95))
    xmax = max(xmax, .10)
    grid = np.linspace(0, xmax, 200)
    logrr = grid / .1 * r.log_rr
    q = stats.t.ppf(.975, int(r.donors)-1)
    sem = grid / .1 * r.se_cluster
    axes[1].plot(grid*100, np.exp(logrr), color='#26364A', lw=2.5)
    axes[1].fill_between(grid*100, np.exp(logrr-q*sem), np.exp(logrr+q*sem), color='#BAC7D7', alpha=.7)
    axes[1].axhline(1, color='gray', ls='--', lw=1)
    axes[1].set_xlabel('Nearby immune cells allocated inside core (%)')
    axes[1].set_ylabel('Adjusted Beta-density ratio\nrelative to 0% inside')
    axes[1].set_title('B  Adjusted primary association', loc='left', fontweight='bold')
    axes[1].text(.04, .96, f'Per 10 percentage points:\nRR {r.rr:.3f} (95% CI {r.lower:.3f}–{r.upper:.3f})\n'
                 f'Donor-cluster t p = {r.p_cluster_t:.3g}', transform=axes[1].transAxes,
                 va='top', bbox=dict(boxstyle='round,pad=.5', fc='white', ec='#CED5DE'))
    # Rug is descriptive exposure support; capped for legibility, not a sample-weight assumption.
    vals = d.loc[d.immune_localization_fraction <= xmax, 'immune_localization_fraction'].to_numpy()
    axes[1].plot(vals*100, np.repeat(.015, len(vals)), '|', color='#778596', alpha=.15,
                 transform=axes[1].get_xaxis_transform())
    fig.suptitle('Immune localization and current Beta-cell density', fontsize=15, fontweight='bold')
    fig.text(.5, .015, 'A: unadjusted donor means. B: log-linear QMLE with donor, size/class, nearby burden, '
             'band geometry and independent regional-field adjustment. Association is not causation.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .07, 1, .94])
    save_figure(fig, 'figure2_primary_association', out)

    main = result.loc[(result.term == 'allocation_10pp') & (result.outcome == 'n_beta_core')].copy()
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5), gridspec_kw={'width_ratios': [1.5, 1]})
    yy = np.arange(len(main))
    axes[0].errorbar(main.rr, yy, xerr=[main.rr-main.lower, main.upper-main.rr], fmt='o', color='#26364A', capsize=3)
    axes[0].axvline(1, color='gray', ls='--')
    axes[0].set_yticks(yy, [f'{r.model}  (D={r.donors}, n={r.objects})' for r in main.itertuples()])
    axes[0].invert_yaxis()
    axes[0].set_xscale('log')
    axes[0].xaxis.set_major_locator(FixedLocator([.8, .9, 1., 1.1, 1.2]))
    axes[0].xaxis.set_major_formatter(FixedFormatter(['0.80', '0.90', '1.00', '1.10', '1.20']))
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].set_xlabel('Beta-density ratio per 10 percentage points inside')
    axes[0].set_title('A  Prespecified sensitivity analyses', loc='left', fontweight='bold')
    yy2 = np.arange(len(lodo))
    for j, r2 in enumerate(lodo.itertuples()):
        axes[1].scatter(r2.rr, j, color=COLORS[r2.donor_status], s=35)
    axes[1].axvline(r.rr, color='#26364A', label='Primary estimate')
    axes[1].axvline(1, color='gray', ls='--')
    axes[1].set_yticks(yy2, ['Omit '+str(x) for x in lodo.omitted_donor])
    axes[1].invert_yaxis()
    axes[1].set_xscale('log')
    axes[1].xaxis.set_major_locator(FixedLocator([.94, .96, .98, 1.]))
    axes[1].xaxis.set_major_formatter(FixedFormatter(['0.94', '0.96', '0.98', '1.00']))
    axes[1].xaxis.set_minor_locator(NullLocator())
    axes[1].set_xlabel('Refitted density ratio')
    axes[1].set_title('B  Leave one donor out', loc='left', fontweight='bold')
    axes[1].legend(fontsize=9)
    fig.suptitle('How dependent is H1 on the analysis choices and individual donors?', fontsize=15, fontweight='bold')
    fig.text(.5, .015, 'Intervals use donor-cluster t uncertainty. Small subgroup estimates are exploratory; '
             'rows are related sensitivity analyses, not independent replications.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .05, 1, .94])
    save_figure(fig, 'figure3_sensitivity', out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--features', type=Path, default=OUT / 'h1_object_features.csv')
    p.add_argument('--output', type=Path, default=OUT)
    p.add_argument('--bootstrap', type=int, default=999)
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.features, dtype={'case_id': str})
    assert df.case_id.nunique() == 20
    if 'core_complete' not in df:
        assert len(df) == 4763
    assert not df.islet_id.duplicated().any()
    d = eligible(df)
    support = support_table(df, d)
    support.to_csv(args.output / 'donor_analysis_support.csv', index=False)
    d.to_csv(args.output / 'primary_analysis_objects.csv', index=False)
    allrows, model_details, failures = [], {}, []
    primaryrows, primarydetails, fit = fit_model(d, 'Primary', wild=True)
    allrows.extend(primaryrows)
    model_details['Primary'] = primarydetails
    print('Primary model:', json.dumps(primaryrows), flush=True)
    specifications = [
        ('Unweighted objects', d, {'balanced': False}),
        ('Islets only', d.loc[d.object_class == 'Islet'], {}),
        ('At least 10 core cells', d.loc[d.n_cells_core >= 10], {}),
        ('Nearby pool ≥3', d.loc[d.n_immune_nearby >= 3], {}),
        ('Nearby pool ≥5', d.loc[d.n_immune_nearby >= 5], {}),
        ('Aab+ and T1D', d.loc[d.donor_status != 'ND'], {}),
        ('T1D only', d.loc[d.donor_status == 'T1D'], {}),
        ('Original acquisitions', d.loc[~d.case_id.isin(ADDITIONS)], {}),
        ('100–250 µm field', eligible(df, field='100_250'), {}),
        ('All positive field areas', eligible(df, min_field_area=np.nextafter(0., 1.)), {}),
        ('Exclusive peri allocation', eligible(df, exclusive=True), {}),
        ('Two-part localization', d, {'exposure': 'two_part'}),
        ('Non-Beta endocrine outcome', d, {'outcome': 'nonbeta_endocrine'}),
        ('Total core cells outcome', d, {'outcome': 'n_cells_core'}),
    ]
    if 'core_tissue_fraction' in d:
        specifications.append(('Core fully in Tissue', d.loc[d.core_tissue_fraction >= .999], {}))
    for name, sub, options in specifications:
        try:
            rows, detail, _ = fit_model(sub, name, **options)
            allrows.extend(rows)
            model_details[name] = detail
        except Exception as exc:
            failures.append(dict(model=name, error=str(exc)))
    result = pd.DataFrame(allrows)
    result.to_csv(args.output / 'model_results.csv', index=False)
    write_json(args.output / 'model_details.json', model_details)
    write_json(args.output / 'model_failures.json', failures)
    lodo = []
    for donor in sorted(d.case_id.unique()):
        try:
            rows, _, _ = fit_model(d.loc[d.case_id != donor], 'LODO')
            lodo.append(dict(omitted_donor=donor, donor_status=d.loc[d.case_id == donor, 'donor_status'].iloc[0], **rows[0]))
        except Exception as exc:
            failures.append(dict(model='LODO_'+donor, error=str(exc)))
    lodo = pd.DataFrame(lodo)
    lodo.to_csv(args.output / 'leave_one_donor_out.csv', index=False)
    figures(df, d, support, result, model_details, fit, lodo, args.output)
    if args.bootstrap:
        boots = []
        boot_data = d[['case_id', 'object_class', 'core_area_um2', 'peri_area_um2',
                       'n_immune_nearby', 'immune_localization_fraction', 'field_density',
                       'n_beta_core', 'has_core_immune']].copy()
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for i, r in enumerate(pool.map(boot_one, [(boot_data, k) for k in range(args.bootstrap)], chunksize=5)):
                boots.append(r)
                if (i+1) % 100 == 0:
                    print(f'Donor bootstrap {i+1}/{args.bootstrap}', flush=True)
        bdf = pd.DataFrame(boots)
        bdf.to_csv(args.output / 'donor_bootstrap.csv', index=False)
        valid = bdf.log_rr.dropna()
        boot_result = dict(requested=args.bootstrap, completed=len(valid), failed=args.bootstrap-len(valid),
                           rr_percentile_ci=np.exp(valid.quantile([.025, .975]).to_numpy()).tolist(),
                           median_rr=float(np.exp(valid.median())), seed=SEED)
    else:
        boot_result = None
    summary = dict(data_root='data/development/uniform-20donor-20261022-v10',
        edge_eligibility_applied='core_complete' in df,
        core_complete_objects=int(df.core_complete.sum()) if 'core_complete' in df else None,
        features_sha256=hashlib.sha256(args.features.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        paired_objects=len(df), donors=int(df.case_id.nunique()), eligible_objects=len(d),
        eligible_core_positive_objects=int(d.has_core_immune.sum()),
        eligible_core_immune_cells=int(d.n_immune_core.sum()),
        no_nearby_immune_objects=int((df.n_immune_nearby == 0).sum()),
        primary=primaryrows[0], wild_score=primarydetails.get('wild_restricted_score'),
        donor_bootstrap=boot_result, model_failures=failures,
        lodo_rr_range=[float(lodo.rr.min()), float(lodo.rr.max())],
        versions={'numpy': np.__version__, 'pandas': pd.__version__, 'statsmodels': sm.__version__},
        interpretation='Cross-sectional association, not killing or cumulative injury; exploratory, donor replication limited')
    write_json(args.output / 'analysis_summary.json', summary)
    write_json(args.output / 'model_failures.json', failures)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
