#!/usr/bin/env python3
"""Analysis-only H1 follow-ups on the frozen corrected object table."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/h1-followup-mpl')
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import h1_analysis_20261022 as h1

OUT = Path('output/h1_followups_h2_h4_20261024')
SOURCE = Path('output/h1_edge_corrected_20261024/h1_object_features.parquet')


def prepare(d, kind):
    d = d.copy()
    dropped = []
    if kind == 'presence':
        d['present'] = (d.n_beta_core > 0).astype(int)
        varied = d.groupby('case_id').present.nunique()
        dropped = varied[varied < 2].index.tolist()
        d = d.loc[~d.case_id.isin(dropped)].copy()
    elif kind == 'positive_density':
        d = d.loc[d.n_beta_core > 0].copy()
    return d, dropped


def fit(d, kind):
    X = h1.build_design(d)
    term = 'allocation_10pp'
    if kind == 'reciprocal':
        X = X.drop(columns=term)
        term = 'beta_density_log2_shifted'
        X[term] = np.log2(1 + d.beta_density_mm2 / 1000)
        y = np.column_stack([d.n_immune_core, d.n_immune_nearby-d.n_immune_core])
        assert (y >= 0).all() and (y.sum(axis=1) > 0).all()
        offset = np.log(d.core_area_um2 / d.peri_area_um2)
        family = sm.families.Binomial()
    elif kind == 'presence':
        y = d.present
        offset = np.zeros(len(d))
        family = sm.families.Binomial()
    elif kind == 'positive_density':
        y = d.beta_density_mm2
        offset = np.zeros(len(d))
        family = sm.families.Gamma(link=sm.families.links.Log())
    elif kind == 'count':
        y = d.n_beta_core
        offset = np.log(d.core_area_um2 / 1e6)
        family = sm.families.Poisson()
    else:
        raise ValueError(kind)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError('Rank deficient design')
    w = 1 / d.case_id.map(d.groupby('case_id').size()).to_numpy(float)
    w /= w.mean()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fitted = sm.GLM(y, X, family=family, offset=offset, freq_weights=w).fit(
            maxiter=500, tol=1e-10)
    messages = sorted(set(str(c.message) for c in caught))
    if not fitted.converged or not np.isfinite(fitted.params).all():
        raise ValueError('Mean model failed to converge')
    if any('separation' in s.lower() for s in messages):
        raise ValueError('Separation: '+str(messages))
    return float(fitted.params[term]), dict(
        columns=X.columns.tolist(), coefficients=fitted.params.to_dict(),
        warnings=messages, converged=True,
        donor_total_weights={g: float(w[d.case_id.eq(g)].sum()) for g in sorted(d.case_id.unique())})


def jackknife(d, kind, group, variant='shared_peri'):
    data, dropped = prepare(d, kind)
    donors = sorted(data.case_id.unique())
    base = dict(group=group, model=kind, variant=variant, donors=len(donors),
                objects=len(data), beta_positive=int((data.n_beta_core > 0).sum()),
                excluded_constant_binary_donors=','.join(dropped))
    if len(donors) < 3:
        return dict(base, status='insufficient donors'), [], {}
    try:
        b, details = fit(data, kind)
    except Exception as e:
        return dict(base, status='fit_failed', error=str(e)), [], {}
    deletes = []
    for g in donors:
        row = dict(group=group, model=kind, variant=variant, deleted_donor=g)
        try:
            value, _ = fit(data.loc[data.case_id.ne(g)].copy(), kind)
            row.update(log_effect=value, status='ok')
        except Exception as e:
            row.update(log_effect=np.nan, status=str(e))
        deletes.append(row)
    estimates = np.array([r['log_effect'] for r in deletes])
    if not np.isfinite(estimates).all():
        return dict(base, status='jackknife_failed', log_effect=b, effect=np.exp(b)), deletes, details
    se = np.sqrt((len(donors)-1)/len(donors) * np.sum((estimates-estimates.mean())**2))
    df = len(donors)-1
    q = stats.t.ppf(.975, df)
    return dict(base, status='ok', log_effect=b, effect=np.exp(b), se_jackknife=se,
                df=df, lower=np.exp(b-q*se), upper=np.exp(b+q*se),
                p=2*stats.t.sf(abs(b/se), df),
                delete_min=np.exp(estimates.min()), delete_max=np.exp(estimates.max())), deletes, details


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_parquet(SOURCE)
    rows, deletes, details = [], [], {}
    for exclusive in [False, True]:
        d = h1.eligible(raw, exclusive=exclusive)
        variant = 'exclusive_peri' if exclusive else 'shared_peri'
        for group in ['Pooled', *h1.STATUS]:
            sub = d if group == 'Pooled' else d.loc[d.donor_status.eq(group)]
            for kind in (['reciprocal'] if exclusive else ['count', 'presence', 'positive_density', 'reciprocal']):
                r, j, detail = jackknife(sub, kind, group, variant)
                rows.append(r); deletes.extend(j)
                details[f'{variant}/{group}/{kind}'] = detail
                print(json.dumps(r), flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(OUT/'h1_followup_models.csv', index=False)
    pd.DataFrame(deletes).to_csv(OUT/'h1_donor_deletion_fits.csv', index=False)
    (OUT/'h1_model_details.json').write_text(json.dumps(details, indent=2))
    contrasts = []
    count = result.loc[result.model.eq('count') & result.group.ne('Pooled')].set_index('group')
    for a,b in [('T1D','ND'), ('T1D','Aab+'), ('Aab+','ND')]:
        ra,rb = count.loc[a],count.loc[b]
        va,vb = ra.se_jackknife**2,rb.se_jackknife**2
        delta = ra.log_effect-rb.log_effect
        se = np.sqrt(va+vb)
        df = (va+vb)**2/(va**2/ra.df+vb**2/rb.df)
        q=stats.t.ppf(.975,df)
        contrasts.append(dict(contrast=f'{a} vs {b}', ratio_of_RRs=np.exp(delta),
                              lower=np.exp(delta-q*se), upper=np.exp(delta+q*se),
                              df=df, p=2*stats.t.sf(abs(delta/se),df)))
    c=pd.DataFrame(contrasts);c['p_holm']=multipletests(c.p,method='holm')[1]
    c.to_csv(OUT/'h1_group_contrasts.csv',index=False)
    d=h1.eligible(raw)
    census=d.assign(beta_present=d.n_beta_core.gt(0)).groupby(
        ['donor_status','case_id','has_core_immune']).agg(
        objects=('islet_id','size'), beta_positive=('beta_present','sum'),
        beta_cells=('n_beta_core','sum')).reset_index()
    census['beta_presence_fraction']=census.beta_positive/census.objects
    census.to_csv(OUT/'h1_beta_presence_census.csv',index=False)
    plot(result)


def plot(result):
    fig,axs=plt.subplots(2,2,figsize=(12,8))
    labels={'count':'Beta count / area: density ratio',
            'presence':'Any Beta cells: odds ratio',
            'positive_density':'Beta-positive sections: density ratio',
            'reciprocal':'Reciprocal: core immune odds ratio'}
    colors=['#555555',*h1.COLORS.values()]
    for ax,kind in zip(axs.flat, labels):
        r=result.loc[result.model.eq(kind)&result.variant.eq('shared_peri')]
        for i,(_,row) in enumerate(r.iterrows()):
            ax.errorbar(row.effect,3-i,xerr=[[row.effect-row.lower],[row.upper-row.effect]],
                        fmt='o',color=colors[i],capsize=4)
            ax.text(.02,3-i-.31,f'{row.effect:.2f} [{row.lower:.2f}, {row.upper:.2f}]',
                    transform=ax.get_yaxis_transform(),fontsize=9,va='center')
        ax.axvline(1,color='#999999',linestyle='--');ax.set_xscale('log')
        ax.set_yticks(range(4),['T1D (6 donors)','Aab+ (7 donors)','ND (7 donors)','Pooled (20 donors)'])
        ax.set_ylim(-.6,3.8);ax.set_title(labels[kind],fontsize=12)
        ticks=[v for v in [.4,.5,.75,1,1.25,1.5,2] if ax.get_xlim()[0]<=v<=ax.get_xlim()[1]]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FixedFormatter([str(v) for v in ticks]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlabel('Per doubling of (1 + Beta density / 1000)' if kind=='reciprocal'
                      else 'Per 10 percentage points of core immune allocation')
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('H1 follow-ups on the same corrected data\nWhole-donor jackknife 95% intervals; exploratory inference',fontsize=15)
    fig.tight_layout(rect=(0,0,1,.93))
    fig.savefig(OUT/'h1_three_checks.png',dpi=180);fig.savefig(OUT/'h1_three_checks.pdf')


if __name__=='__main__':
    run()
