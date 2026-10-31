#!/usr/bin/env python3
"""Supporting original-direction audit of existing independent donor coefficients.

No cell/object model is refitted. This audit supplements the preserved primary
disease comparisons. Branch decisions are unassigned (user judgment).
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h1-h2-h5-conditional-audit-20261026')
from pathlib import Path
import hashlib,json,re,sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h1_h10_conditional_audit_20261026/association_audit'
GROUPS=['ND','Aab+','T1D']
COLORS=['#4477AA','#CC6633','#228833']
NARRATIVE_ONLY='--narrative-only' in sys.argv
SPECS={
 'H1':dict(source=ROOT/'output/h1_h4_disease_comparisons_20261026/h1/donor_coefficients.csv',
  title='Immune localization\nand β-cell density',
  exposure='Fraction of nearby immune cells located within the endocrine core',
  outcome='Core Beta-cell density',
  controls='log core area, log peri/core area ratio, log1p nearby immune count, log1p regional immune density, and object class'),
 'H2':dict(source=ROOT/'output/h1_h4_disease_comparisons_20261026/h2/primary_donor_coefficients.csv',
  title='ColIV coverage\nand core immune density',
  exposure='Nonvascular ColIV-positive coverage in the 2-µm boundary band',
  outcome='Core immune-cell density',
  controls='log1p peri immune density, log core area, log perimeter, log1p tissue-edge distance, capillary area fraction, and object class'),
 'H5':dict(source=ROOT/'output/h5_h10_original_plan_20261026/h5/plot_data_donor_coefficients.csv',
  title='Immune localization\nand β-cell INS',
  exposure='Fraction of nearby immune cells located within the endocrine core',
  outcome='Mean scaled INS within retained typed Beta cells',
  controls='log1p Beta-cell density, log core area, object class; donor-constant background absorbed by the donor intercept')}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(values,hypothesis,group):
    x=np.asarray(values,float);assert len(x)>=2 and np.isfinite(x).all()
    n=len(x);mean=x.mean();sd=x.std(ddof=1);se=sd/np.sqrt(n)
    test=stats.ttest_1samp(x,popmean=0,alternative='two-sided')
    delta=stats.t.ppf(.975,n-1)*se
    return dict(hypothesis=hypothesis,donor_status=group,donors=n,mean_beta=mean,sd_between_donors=sd,
        SE=se,lower_95=mean-delta,upper_95=mean+delta,df=n-1,t_statistic=test.statistic,
        p_value=test.pvalue,negative_donors=int((x<0).sum()),positive_donors=int((x>0).sum()),
        zero_donors=int((x==0).sum()),expected_direction='Negative',
        analysis_role='Supporting original-direction audit added on2026-10-26; not a new primary',
        test='Two-sided one-sample t test of mean independent donor coefficient versus zero')


def verify_h10():
    h10=ROOT/'output/h5_h10_original_plan_20261026/h10'
    h5=ROOT/'output/h5_h10_original_plan_20261026/h5/plot_data_donor_ins_measurements.csv'
    donor=pd.read_csv(h10/'donor_inputs.csv').set_index('case_id')
    measured=donor[~donor.C_peptide_censored]
    assert len(donor)==20 and donor.index.is_unique and len(measured)==19 and 78 not in donor.index
    assert donor.loc[donor.C_peptide_censored].index.tolist()==[6550]
    checks=[]
    for prefix,expected in [('clinical',measured),('clinical_disease_conditioned',measured)]:
        a=pd.read_csv(h10/f'{prefix}_held_out_predictions.csv').set_index('case_id')
        assert a.index.is_unique and set(a.index)==set(expected.index)
        assert np.allclose(a.observed,np.log(expected.loc[a.index,'C_peptide_ng_mL']),rtol=1e-14,atol=1e-14)
        assert a[['baseline_prediction','augmented_prediction']].notna().all().all()
        for model in ['baseline','augmented']:
            assert np.allclose(a[model+'_squared_error'],(a[model+'_prediction']-a.observed)**2)
        checks.append(dict(export=prefix+'_held_out_predictions.csv',rows=len(a),outcome='Natural log measured C-peptide',complete=True))
    source=pd.read_csv(h5).set_index('case_id')
    inp=pd.read_csv(h10/'tissue_followup_inputs.csv').set_index('case_id')
    assert set(inp.index)==set(source.index)==set(donor.index)
    assert np.allclose(inp.mean_beta_INS,source.loc[inp.index,'mean_beta_INS'],rtol=1e-14,atol=1e-14)
    assert np.array_equal(inp.n_beta_ins_measured,source.loc[inp.index,'n_beta_ins_measured'])
    assert inp.n_beta_ins_measured.gt(0).all()
    for prefix in ['tissue_INS','tissue_INS_disease_conditioned']:
        a=pd.read_csv(h10/f'{prefix}_held_out_predictions.csv').set_index('case_id')
        assert a.index.is_unique and set(a.index)==set(source.index)
        assert np.allclose(a.observed,source.loc[a.index,'mean_beta_INS'],rtol=1e-14,atol=1e-14)
        assert a[['baseline_prediction','augmented_prediction']].notna().all().all()
        checks.append(dict(export=prefix+'_held_out_predictions.csv',rows=len(a),outcome='Existing H5 cell-weighted mean INS within typed Beta cells',complete=True))
    t=pd.read_csv(h10/'all_six_T1D_observations_and_predictions.csv').set_index('case_id')
    assert len(t)==6 and set(t.index)==set(donor.index[donor.donor_status.eq('T1D')])
    assert t.loc[6550,'C_peptide_censored'] and t.loc[6550,'C_peptide_ng_mL']==.02
    assert pd.isna(t.loc[6550,'baseline_log_error']) and pd.isna(t.loc[6550,'augmented_log_error'])
    assert t[['baseline_prediction','augmented_prediction']].notna().all().all()
    cens=pd.read_csv(h10/'censored_model_fitted_values.csv')
    for _,s in cens.groupby(['analysis','model']):
        assert len(s)==20 and s.case_id.nunique()==20 and set(s.case_id)==set(donor.index)
        assert s.C_peptide_censored.sum()==1 and s.loc[s.C_peptide_censored,'case_id'].tolist()==[6550]
    plots={'clinical_held_out_prediction':19,'clinical_disease_conditioned_prediction':19,
        'tissue_INS_held_out_prediction':20,'tissue_INS_disease_conditioned_prediction':20,
        'T1D_individual_clinical_predictions':6,'censored_model_score_association':2}
    for stem,n in plots.items():
        for suffix in ['.png','.svg','_legend.txt']:assert (h10/(stem+suffix)).is_file()
        frame=pd.read_csv(h10/f'prism_{stem}.csv');assert len(frame)==n
    report=dict(status='Passed: existing outcomes and follow-up exports remain complete; no models refitted',
        primary_measured_donors=19,clinical_outcome='ln measured C-peptide',censored_donors=1,censored_case=6550,
        all_T1D_donors_displayed=6,tissue_followup_donors=20,
        tissue_followup_outcome='Exact existing H5 cell-weighted mean Beta-cell INS; no invented zeros',
        plot_csv_rows=plots,prediction_checks=checks,
        censored_fits_preserve20_donors=True,censored_donor_squared_error_is_missing=True,
        source_H5_sha256=sha(h5),source_H10_inputs_sha256=sha(h10/'donor_inputs.csv'))
    (OUT/'H10_export_completeness.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    summaries=[];all_summaries=[];donor_rows=[];provenance={}
    for hypothesis,spec in SPECS.items():
        path=spec['source'];d=pd.read_csv(path)
        assert d.case_id.is_unique and 78 not in d.case_id.values
        d=d[np.isfinite(d.standardized_beta)].copy()
        d['hypothesis']=hypothesis;donor_rows.append(d)
        assert d.donor_status.isin(GROUPS).all()
        if NARRATIVE_ONLY:
            existing=pd.read_csv(OUT/'state_conditional_means.csv')
            rows=existing[existing.hypothesis.eq(hypothesis)].to_dict('records')
        else:
            rows=[summarize(d.loc[d.donor_status.eq(group),'standardized_beta'],hypothesis,group) for group in GROUPS]
            q=multipletests([r['p_value'] for r in rows],method='fdr_bh')[1]
            for r,p in zip(rows,q):r['p_bh_three_states']=p
        summaries+=rows
        if NARRATIVE_ONLY:
            existing=pd.read_csv(OUT/'equal_donor_all_donor_means.csv')
            pooled=existing[existing.hypothesis.eq(hypothesis)].iloc[0].to_dict()
        else:pooled=summarize(d.standardized_beta,hypothesis,'All donors equally weighted')
        pooled['p_bh_three_states']=np.nan
        pooled['analysis_role']='Supporting all-donor mean only; each donor equally weighted; not a disease contrast or primary test'
        all_summaries.append(pooled)
        provenance[hypothesis]=dict(source=str(path),source_sha256=sha(path),donors=len(d),
            exposure=spec['exposure'],outcome=spec['outcome'],controls=spec['controls'],cell_models_refitted=False)
    summary=pd.DataFrame(summaries);pooled=pd.DataFrame(all_summaries);donors=pd.concat(donor_rows,ignore_index=True)
    if not NARRATIVE_ONLY:
        summary.to_csv(OUT/'state_conditional_means.csv',index=False)
        pooled.to_csv(OUT/'equal_donor_all_donor_means.csv',index=False)
        donors.to_csv(OUT/'exact_donor_coefficients.csv',index=False)
        summary.to_csv(OUT/'plot_data_conditional_means.csv',index=False)
    plots=[]
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,3,figsize=(12.5,4.5))
    for ax,(hypothesis,spec) in zip(axs,SPECS.items()):
        rows=summary[summary.hypothesis.eq(hypothesis)].set_index('donor_status').loc[GROUPS]
        lo=min(rows.lower_95.min(),0);hi=max(rows.upper_95.max(),0);span=max(hi-lo,.01)
        for i,(group,color) in enumerate(zip(GROUPS,COLORS)):
            r=rows.loc[group]
            ax.errorbar(i,r.mean_beta,yerr=[[r.mean_beta-r.lower_95],[r.upper_95-r.mean_beta]],
                fmt='o',color=color,ms=7,lw=1.5,capsize=4)
            if r.p_bh_three_states<.05:
                star='***' if r.p_bh_three_states<.001 else '**' if r.p_bh_three_states<.01 else '*'
                ax.text(i,r.upper_95+.035*span,star,ha='center',va='bottom',fontsize=14)
        ax.axhline(0,color='#888888',ls='--',lw=.9)
        ax.set_xticks(range(3),GROUPS);ax.set_xlim(-.5,2.5);ax.set_ylim(lo-.08*span,hi+.17*span)
        ax.set_title(spec['title'],fontsize=11);ax.set_xlabel('Disease state')
        ax.spines[['top','right']].set_visible(False)
    axs[0].set_ylabel('Mean adjusted donor coefficient (β)\n95% confidence interval')
    if not NARRATIVE_ONLY:
        fig.tight_layout();fig.savefig(OUT/'conditional_association_means.png',dpi=220)
        fig.savefig(OUT/'conditional_association_means.svg')
    plt.close(fig)
    text=('Supporting audit of the original directional associations. Points show the arithmetic mean of independent adjusted donor regression coefficients; '
        'bars show pointwise 95% Student t confidence intervals using n−1 degrees of freedom. Every donor receives equal weight, irrespective of its number of objects. '
        'Left: immune localization versus Beta-cell density; middle: nonvascular ColIV boundary coverage versus core immune density; right: immune localization versus '
        'mean INS in retained typed Beta cells. Exposure and continuous covariates were standardized within donor, and each outcome across its included objects, '
        'as in the existing source models. These coefficients are not correlations and the three panels represent different outcomes. '
        'ND and Aab+ each contribute seven donors in every panel; T1D contributes six in the first two panels and five in the INS panel. '
        'Donor6563 has no estimable INS association in the existing H5 donor model and is not assigned zero. '
        'Two-sided one-sample t tests compare each state mean coefficient with zero. Benjamini–Hochberg correction covers the three state-specific tests separately '
        'within each hypothesis. Asterisks mark slopes differing from zero after that correction (*p<0.05, **p<0.01, ***p<0.001); '
        'they do not denote differences between disease states. The displayed 95% intervals are pointwise, not multiplicity-adjusted. '
        'The original predicted direction is negative in all three panels. This supporting analysis was added for the current if–then audit; '
        'it does not replace the preserved primary disease comparisons. Branch decisions are unassigned (user judgment). '
        'The separate all-donor summaries are supporting equal-donor means only and are not plotted.\n\n')
    for hypothesis,spec in SPECS.items():
        text+=f'{spec["title"].replace(chr(10)," ")}: exposure={spec["exposure"]}; outcome={spec["outcome"]}; controls={spec["controls"]}.\n'
        for r in summary[summary.hypothesis.eq(hypothesis)].itertuples():
            text+=f'{r.donor_status}: n={r.donors}, mean β={r.mean_beta:+.5f}, 95% CI [{r.lower_95:+.5f}, {r.upper_95:+.5f}], '
            text+=f't({r.df})={r.t_statistic:+.4f}, unadjusted p={r.p_value:.6g}, BH p={r.p_bh_three_states:.6g}; {r.negative_donors} negative/{r.positive_donors} positive/{r.zero_donors} zero coefficients.\n'
        text+='\n'
    (OUT/'figure_legend.txt').write_text(text)
    report=['Supporting original-direction audit of existing donor coefficients','',
        'No cell/object model was refitted. Existing disease comparisons remain unchanged. All tests here were added for the current if–then audit, not declared as new primary tests.','',
        'State-specific mean adjusted coefficients (pointwise95%t intervals; BH across three states within each hypothesis):']
    for r in summary.itertuples():report.append(f'{r.hypothesis} {r.donor_status}: n={r.donors}; mean={r.mean_beta:+.6f} [{r.lower_95:+.6f}, {r.upper_95:+.6f}]; p={r.p_value:.6f}; BH={r.p_bh_three_states:.6f}; negative/positive/zero={r.negative_donors}/{r.positive_donors}/{r.zero_donors}.')
    report+=['','Equal-donor all-donor means (supporting only; unadjusted p):']
    for r in pooled.itertuples():report.append(f'{r.hypothesis}: n={r.donors}; mean={r.mean_beta:+.6f} [{r.lower_95:+.6f}, {r.upper_95:+.6f}]; p={r.p_value:.6f}.')
    report+=['','Methods: these summaries assume approximately normal independent donor coefficients, '
        'and equal donor weighting differs from pooling all objects or inverse-variance weighting. H5 is conditional on the five estimable T1D donors.','',
        'Original if–then conditions are preserved verbatim in original_if_then_conditions.json. Branch decisions in if_then_decisions.json are unassigned (user judgment).','',
        'Method precedent: repository literature/Heide.pdf, PDF p32, one-sample t tests against a hypothetical mean; this is the same conventional donor-coefficient mean/interval approach used in H7. '
        'The precise supporting role and existing exposure/outcome adjustments are recorded in the separate legend.']
    (OUT/'RESULTS.txt').write_text('\n'.join(report)+'\n')
    plan_path=ROOT/'docs/T1D_research_hypotheses_20261022.md';plan=plan_path.read_text();conditions=[]
    for hypothesis in SPECS:
        section=plan.split('### '+hypothesis+' ',1)[1].split('\n### ',1)[0]
        entries=[dict(condition=title,text=text.strip()) for title,text in re.findall(r'\*\*(If[^*]+)\*\*\s*(.*?)(?=\n\n|\Z)',section,flags=re.S)]
        conditions.append(dict(hypothesis=hypothesis,original_conditions=entries))
    (OUT/'original_if_then_conditions.json').write_text(json.dumps(dict(source=str(plan_path),hypotheses=conditions),indent=2)+'\n')
    (OUT/'if_then_decisions.json').write_text(json.dumps([dict(hypothesis=h,branch_decision='unassigned (user judgment)') for h in SPECS],indent=2)+'\n')
    (OUT/'provenance.json').write_text(json.dumps(dict(analysis='Supporting current if–then audit',sources=provenance,
        method='One-sample t mean/95%CI; two-sided; BH across3states within each hypothesis',
        original_disease_comparisons_preserved=True,cell_models_refitted=False,script_sha256=sha(__file__)),indent=2)+'\n')
    h10=json.loads((OUT/'H10_export_completeness.json').read_text()) if NARRATIVE_ONLY else verify_h10()
    print('\n'.join(report));print('\nH10 export verification:',h10['status'])


if __name__=='__main__':main()
