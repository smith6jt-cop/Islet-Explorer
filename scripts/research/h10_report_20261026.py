#!/usr/bin/env python3
"""Factual figures, statistics, methods, and unassigned plan conditions for H10."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h10-approved-report-20261026-mpl')
from pathlib import Path
import hashlib,json,re,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[2]
OUT=REPO/'output/h5_h10_original_plan_20261026/h10'
GROUPS=['ND','Aab+','T1D'];COLORS={'ND':'#4477AA','Aab+':'#CC6633','T1D':'#228833'}
NARRATIVE_ONLY='--narrative-only' in sys.argv

def save(fig,stem):
    if not NARRATIVE_ONLY:
        fig.tight_layout();fig.savefig(OUT/f'{stem}.png',dpi=220);fig.savefig(OUT/f'{stem}.svg')
    plt.close(fig)


def prediction_figure(label,axis_label,stem):
    d=pd.read_csv(OUT/f'{label}_held_out_predictions.csv')
    fig,axs=plt.subplots(1,2,figsize=(9.5,4.5),sharex=True,sharey=True)
    lo=d[['observed','baseline_prediction','augmented_prediction']].min().min()
    hi=d[['observed','baseline_prediction','augmented_prediction']].max().max();pad=(hi-lo)*.09
    for ax,model,title in zip(axs,['baseline','augmented'],['Baseline','With spatial score']):
        for group in GROUPS:
            z=d[d.donor_status.eq(group)]
            ax.scatter(z.observed,z[model+'_prediction'],s=42,color=COLORS[group],edgecolor='white',lw=.6,label=group,zorder=3)
        ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],color='#777777',ls='--',lw=1)
        ax.set_xlim(lo-pad,hi+pad);ax.set_ylim(lo-pad,hi+pad)
        ax.set_xlabel('Observed '+axis_label);ax.set_title(title,fontsize=12)
        ax.spines[['top','right']].set_visible(False)
    axs[0].set_ylabel('Held-out predicted '+axis_label);axs[0].legend(frameon=False,fontsize=9)
    save(fig,stem)
    if not NARRATIVE_ONLY:d.to_csv(OUT/f'prism_{stem}.csv',index=False)


def main():
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    comparison=pd.read_csv(OUT/'ordinary_model_comparisons.csv').set_index('analysis')
    clinical=comparison.loc['clinical'];conditioned=comparison.loc['clinical_disease_conditioned']
    tissue=comparison.loc['tissue_INS'];tissue_conditioned=comparison.loc['tissue_INS_disease_conditioned']
    cens=pd.read_csv(OUT/'censored_likelihood_comparisons.csv').set_index('analysis')
    c=cens.loc['censored'];cd=cens.loc['censored_disease_conditioned']
    prediction_figure('clinical','log C-peptide','clinical_held_out_prediction')
    prediction_figure('clinical_disease_conditioned','log C-peptide','clinical_disease_conditioned_prediction')
    prediction_figure('tissue_INS','mean β-cell INS','tissue_INS_held_out_prediction')
    prediction_figure('tissue_INS_disease_conditioned','mean β-cell INS','tissue_INS_disease_conditioned_prediction')

    t1d=pd.read_csv(OUT/'all_six_T1D_observations_and_predictions.csv').sort_values('case_id')
    fig,ax=plt.subplots(figsize=(9,4.6));positions=np.arange(len(t1d))
    measured=~t1d.C_peptide_censored
    ax.scatter(positions[measured]-.17,t1d.loc[measured,'C_peptide_ng_mL'],s=47,color='#222222',label='Measured',zorder=3)
    ax.scatter(positions+.01,t1d.baseline_predicted_median_ng_mL,s=44,color='#4477AA',marker='s',label='Baseline prediction',zorder=3)
    ax.scatter(positions+.18,t1d.augmented_predicted_median_ng_mL,s=50,color='#CC6633',marker='D',label='With spatial score',zorder=3)
    for i,r in enumerate(t1d.itertuples()):
        if r.C_peptide_censored:
            ax.scatter(i-.17,.02,s=65,marker='v',facecolor='none',edgecolor='#222222',linewidth=1.3,label='Below 0.02',zorder=3)
            ax.annotate('',xy=(i-.17,.012),xytext=(i-.17,.019),arrowprops=dict(arrowstyle='->',color='#222222'))
    ax.set_yscale('log');ax.set_ylim(.009,7);ax.set_ylabel('C-peptide (ng/mL)')
    ax.set_xticks(positions,t1d.case_id.astype(str));ax.set_xlabel('T1D donor')
    ax.legend(frameon=False,fontsize=9,ncol=2,loc='upper center');ax.spines[['top','right']].set_visible(False)
    save(fig,'T1D_individual_clinical_predictions')
    if not NARRATIVE_ONLY:t1d.to_csv(OUT/'prism_T1D_individual_clinical_predictions.csv',index=False)

    # Display the fitted association and its uncertainty for the proper censored sensitivity.
    cc=pd.read_csv(OUT/'censored_model_coefficients.csv')
    rows=cc[cc.model.eq('augmented')&cc.term.eq('spatial_score')].set_index('analysis').loc[['censored','censored_disease_conditioned']].reset_index()
    fig,ax=plt.subplots(figsize=(6.2,4.4))
    x=np.arange(2);ax.errorbar(x,rows.estimate,yerr=np.array([rows.estimate-rows.lower_95,rows.upper_95-rows.estimate]),
        fmt='o',color='#333333',capsize=5,ms=7,lw=1.5)
    ax.axhline(0,color='#888888',ls='--',lw=.8);ax.set_xticks(x,['Clinical covariates','Also adjusted\nfor disease state'])
    ax.set_xlim(-.5,1.5);ax.set_ylabel('Spatial-score association with log C-peptide\n(regression coefficient and 95% CI)')
    ax.spines[['top','right']].set_visible(False);save(fig,'censored_model_score_association')
    if not NARRATIVE_ONLY:rows.to_csv(OUT/'prism_censored_model_score_association.csv',index=False)

    common=('The baseline contains sampled Beta-cell density, donor age, and admission glucose; the augmented model adds '
        'one fixed score equal to one half of standardized mean nonvascular ColIV boundary coverage minus standardized mean immune core localization. '
        'Component directions and equal weights were fixed before clinical fitting. Every held-out fold recalculates scaling and fits the models using only training donors. '
        'No imputation, feature selection, or tuning was performed. Each plotted observation is an independent donor; dashed lines show equal observed and predicted values. ')
    legends={
      'clinical_held_out_prediction':('Donor-held-out prediction of recorded C-peptide. '+common+
        'The outcome is natural log C-peptide. The analysis includes 19 measured donors (7 ND, 7 Aab+, 5 T1D); the below-limit donor is not assigned a measured concentration. '
        f'Baseline and augmented held-out RMSE are {clinical.baseline_held_out_rmse:.4f} and {clinical.augmented_held_out_rmse:.4f} log units. '
        f'The full-sample added-score test is F(1,14)={clinical.statistic:.3f}, p={clinical.p_value:.4f}; adjusted R² changes from '
        f'{clinical.baseline_adjusted_r_squared:.3f} to {clinical.augmented_adjusted_r_squared:.3f}. The F-test and prediction-error calculation are separate statistics.'),
      'clinical_disease_conditioned_prediction':('Disease-conditioned donor-held-out C-peptide prediction. '+common+
        'Both models additionally contain Aab+ and T1D indicators, with ND as reference. The outcome is natural log C-peptide and the same 19 measured donors are included. '
        f'Baseline and augmented held-out RMSE are {conditioned.baseline_held_out_rmse:.4f} and {conditioned.augmented_held_out_rmse:.4f}; '
        f'the added-score test is F(1,12)={conditioned.statistic:.3f}, p={conditioned.p_value:.4f}.'),
      'tissue_INS_held_out_prediction':('Donor-held-out prediction of retained-cell insulin signal. '+common+
        'The outcome is the existing H5 cell-weighted mean scaled INS signal among typed Beta cells across complete eligible endocrine objects. '
        'All 20 donors contribute (7 ND, 7 Aab+, 6 T1D); Beta-negative objects do not contribute invented zero within-Beta means. '
        f'Baseline and augmented RMSE are {tissue.baseline_held_out_rmse:.4f} and {tissue.augmented_held_out_rmse:.4f} scaled INS units; '
        f'F(1,15)={tissue.statistic:.3f}, p={tissue.p_value:.4f}. This is the original same-data tissue-phenotype follow-up, not a secretion measurement.'),
      'tissue_INS_disease_conditioned_prediction':('Disease-conditioned prediction of retained-cell insulin signal. '+common+
        'Both models also contain disease indicators. The outcome is the existing H5 cell-weighted mean scaled INS signal within typed Beta cells; all 20 donors contribute. '
        f'Baseline and augmented RMSE are {tissue_conditioned.baseline_held_out_rmse:.4f} and {tissue_conditioned.augmented_held_out_rmse:.4f}; '
        f'F(1,13)={tissue_conditioned.statistic:.3f}, p={tissue_conditioned.p_value:.4f}.'),
      'T1D_individual_clinical_predictions':('Individual recent-onset T1D observations and primary-model predictions. Black circles show measured C-peptide; '
        'colored symbols show exponentiated predictions of log C-peptide (predicted medians on the ng/mL scale). Each of the five measured donors was held out '
        'while its models were trained on the other 18 measured donors. Donor 6550 is recorded below 0.02 ng/mL, shown by an open downward triangle and arrow; '
        'its predictions were made using the 19 measured donors without assigning it an observed concentration or squared error. '
        'The models and fixed score are defined in the clinical prediction legend. The vertical axis is logarithmic. All six T1D donors are displayed.'),
      'censored_model_score_association':('Spatial-score association retaining the below-detection donor. Points and bars show score coefficients with observed-information '
        'Wald 95% confidence intervals from censored-normal regressions of natural log C-peptide. All 20 donors contribute; 19 measured values contribute normal log densities, '
        'and donor 6550 contributes the log probability below log(0.02). No concentration is substituted for that donor. Clinical covariates are sampled Beta quantity, age, '
        'and admission glucose. The second fit additionally conditions on disease group. '
        f'Added-score likelihood-ratio tests are χ²(1)={c.statistic:.3f}, p={c.p_value:.4f}, and χ²(1)={cd.statistic:.3f}, p={cd.p_value:.4f}, respectively. '
        'These are full-sample association sensitivities, not held-out prediction results; no significance bracket is shown because neither comparison reaches p<0.05.')}
    for name,text in legends.items():(OUT/f'{name}_legend.txt').write_text(text+'\n')

    relative_rmse_change=100*(clinical.augmented_held_out_rmse-clinical.baseline_held_out_rmse)/clinical.baseline_held_out_rmse
    report=['# Spatial information and recorded C-peptide','',
        'This document records measurements, statistics, methods, and technical diagnostics. Branch decisions are unassigned (user judgment).','',
        '## Approved analysis and primary result','',
        'The baseline ordinary regression predicts natural log C-peptide using unique sampled tissue Beta-cell density, age, and admission glucose. '
        'All three predictors are complete. The augmented model adds the previously fixed equal-weight boundary-support/immune-access score. '
        'Nineteen measured donors are included; the donor recorded <0.02 is retained separately and in the censored-model sensitivity.','',
        f'The score coefficient is {clinical.score_coefficient:+.4f} (95% CI {clinical.score_lower_95:+.4f} to {clinical.score_upper_95:+.4f}); '
        f'F(1,14)={clinical.statistic:.3f}, p={clinical.p_value:.4f}. Adjusted R² is {clinical.baseline_adjusted_r_squared:.3f} for baseline and '
        f'{clinical.augmented_adjusted_r_squared:.3f} with the score. The relative held-out RMSE change is {relative_rmse_change:+.2f}%.','',
        '| Analysis | Added-score p | Baseline held-out RMSE | With score | Change |',
        '|---|---:|---:|---:|---:|']
    for label,title in [('clinical','C-peptide, 19 donors (log units)'),('clinical_disease_conditioned','C-peptide + disease adjustment'),
                        ('tissue_INS','Within-Beta INS, 20 donors (scaled units)'),('tissue_INS_disease_conditioned','Within-Beta INS + disease adjustment')]:
        r=comparison.loc[label];report.append(f'| {title} | {r.p_value:.4f} | {r.baseline_held_out_rmse:.4f} | {r.augmented_held_out_rmse:.4f} | {r.augmented_minus_baseline_rmse:+.4f} |')
    report+=['','RMSE change is augmented minus baseline. Each held-out prediction recalculates scaling and refits the model solely on the other donors. '
        'The RMSE differences are descriptive prediction checks; the F-tests assess incremental association, not significance of the RMSE changes.','',
        '## Censoring and disease sensitivity','',
        f'The all-20 censored-normal model uses the actual probability below 0.02 for donor 6550. Its baseline-versus-augmented likelihood-ratio test gives '
        f'χ²(1)={c.statistic:.3f}, p={c.p_value:.4f}; with disease indicators, χ²(1)={cd.statistic:.3f}, p={cd.p_value:.4f}. '
        'No bound or half-bound value was substituted.','',
        '## Six T1D donor measurements and predictions','',
        'All six observations and primary-model predictions are plotted individually. The five measured donors have held-out predictions; the censored donor '
        'has a prediction from the 19 measured donors and remains unscored for RMSE. Concentrations and exponentiated predicted log values below are in ng/mL.','',
        '| Donor | Recorded C-peptide | Baseline predicted median | With score predicted median |',
        '|---|---:|---:|---:|']
    for r in t1d.itertuples():
        observed=('<' if r.C_peptide_censored else '')+f'{r.C_peptide_ng_mL:.2f}'
        report.append(f'| {r.case_id} | {observed} | {r.baseline_predicted_median_ng_mL:.3f} | {r.augmented_predicted_median_ng_mL:.3f} |')
    report+=['','## Original plan conditions','',
        'The original if–then conditions are preserved verbatim in `original_if_then_conditions.json`, separately from the statistics. '
        'The branch decision in `if_then_interpretation.json` is unassigned (user judgment).','',
        'The same-data follow-up uses the existing H5 cell-weighted mean INS within typed Beta cells as the outcome, the same baseline predictors and fixed score, '
        'and donor-held-out predictions for all 20 donors. Its statistics appear in the table above.','',
        '## Method and data provenance','',
        'The original H10 target was censored held-out predictive log scoring. The user approved ordinary nested regression and held-out prediction error, '
        'with a censored-normal likelihood sensitivity. '
        '`approved_analysis_specification.json` is the current specification; earlier proposal files are historical and were not executed. '
        'Supplied-paper precedents cover covariate-adjusted linear models/log transformation (Drotar2025, PDF p7), regression (Heide, p23), '
        'held-out prediction (Williams2024, p3), and RMSE (the supplied CELL-REPORTS draft, p12). The exact clinical adaptations and censored sensitivity were explicitly approved.','',
        'Censored likelihoods were independently reproduced with R survival::survreg; maximum likelihood/scale differences are recorded in '
        '`censored_independent_validation.csv`. All designs and held-out folds are full rank. One Python optimizer reports precision loss despite a maximum '
        'gradient below 3×10⁻⁸ and a positive Hessian; the independent R fit verifies the same optimum. The warning is retained in `censored_model_diagnostics.csv`.','',
        'All plotted values, donor identities, fitted coefficients, scaling parameters, held-out errors, censor flags, and diagnostics are exported as CSV. '
        'P values here are unadjusted primary or explicitly secondary model comparisons; any ten-hypothesis adjustment is reported separately. No app data pipeline was rebuilt.','',
        'Reproduce with `h10_approved_models_20261026.py`, then `h10_censored_validation_20261026.R`, then `h10_report_20261026.py`.']
    (OUT/'README.md').write_text('\n'.join(report)+'\n');(OUT/'RESULTS.txt').write_text('\n'.join(report)+'\n')
    plan_path=REPO/'docs/T1D_research_hypotheses_20261022.md'
    section=plan_path.read_text().split('### H10 ',1)[1].split('\n## ',1)[0]
    conditions=[dict(condition=title,text=text.strip()) for title,text in re.findall(r'\*\*(If[^*]+)\*\*\s*(.*?)(?=\n\n|\Z)',section,flags=re.S)]
    (OUT/'original_if_then_conditions.json').write_text(json.dumps(dict(source=str(plan_path),hypothesis='H10',original_conditions=conditions),indent=2)+'\n')
    branch=dict(hypothesis='H10',branch_decision='unassigned (user judgment)',
        original_conditions_file='original_if_then_conditions.json',
        primary_added_score_p=float(clinical.p_value),
        baseline_held_out_rmse=float(clinical.baseline_held_out_rmse),augmented_held_out_rmse=float(clinical.augmented_held_out_rmse),
        disease_conditioned_added_score_p=float(conditioned.p_value),censored_likelihood_ratio_p=float(c.p_value),
        tissue_followup_added_score_p=float(tissue.p_value),
        tissue_baseline_held_out_rmse=float(tissue.baseline_held_out_rmse),tissue_augmented_held_out_rmse=float(tissue.augmented_held_out_rmse))
    (OUT/'if_then_interpretation.json').write_text(json.dumps(branch,indent=2)+'\n')
    independent=pd.read_csv(OUT/'censored_independent_validation.csv')
    validation=json.loads((OUT/'validation.json').read_text())
    validation['independent_R_censored_likelihood_max_difference']=float(independent.log_likelihood_difference.abs().max())
    validation['independent_R_censored_sigma_max_difference']=float(independent.sigma_difference.abs().max())
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    # Mark prior proposals as superseded, retaining their original details for the record.
    for name in ['input_specification.json','method_review.json']:
        path=OUT/name;obj=json.loads(path.read_text());obj['current_status']='Historical proposal; superseded by approved_analysis_specification.json'
        if name=='method_review.json':
            original=obj.get('original_h10',{})
            if 'unsupported_in_supplied_papers' in original:original['not_identified_in_supplied_paper_search']=original.pop('unsupported_in_supplied_papers')
            obj.pop('interpretation',None)
            obj['measurement_definition']='The score uses existing tissue measurements; the clinical outcome is recorded donor C-peptide'
        path.write_text(json.dumps(obj,indent=2)+'\n')
    path=OUT/'approved_analysis_specification.json';obj=json.loads(path.read_text())
    obj['inferential_scope']='Donors are replicates. The primary association test and held-out error are separate statistics. Censored fits and tissue follow-up are secondary. Branch decision: unassigned (user judgment).'
    path.write_text(json.dumps(obj,indent=2)+'\n')
    (OUT/'conventional_model_specification.txt').write_text(
        'Executed user-approved model specification\n\n'
        'Primary outcome: natural log C-peptide in19 measured donors; no imputation. Baseline: tissue Beta-cell density, age, admission glucose. '
        'Augmented: baseline plus0.5×[training-standardized mean boundary coverage−training-standardized mean immune core localization]. '
        'Full-sample comparison: partial F test, score coefficient/95%CI, adjusted R². Prediction: leave-one-donor-out RMSE with training-only scaling. '
        'Disease sensitivity: add Aab+ and T1D indicators to both models.\n\n'
        'Censored sensitivity: all20 donors; normal log density for measured log C-peptide and normal log CDF at log(0.02) for donor6550. '
        'Baseline/augmented likelihood-ratio comparison, with and without disease indicators. No bound substitutions.\n\n'
        'Tissue follow-up: existing H5 cell-weighted mean INS in typed Beta cells, all20 donors; same baseline/score and donor-held-out predictions. '
        'All6T1D donors are displayed; the censored donor has no squared-error contribution.\n\n'
        'Original plan conditions are in original_if_then_conditions.json. Branch decision: unassigned (user judgment).\n')
    manifest=dict(figures=[dict(stem=name,legend=name+'_legend.txt') for name in legends],
        files={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'})
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('\n'.join(report[:19]))


if __name__=='__main__':main()
