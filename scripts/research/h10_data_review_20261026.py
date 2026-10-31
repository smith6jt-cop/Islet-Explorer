#!/usr/bin/env python3
"""Describe H10 measurements and approved-method choices; no clinical model fit."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h10-data-review-20261026-mpl')
from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[2]
OUT=REPO/'output/h5_h10_original_plan_20261026/h10'

def main():
    d=pd.read_csv(OUT/'donor_inputs.csv')
    rows=[]
    for col in ['HbA1c_percent','admission_glucose_mg_dL','peak_glucose_mg_dL']:
        for group,z in [('All',d)]+list(d.groupby('donor_status')):
            rows.append(dict(variable=col,group=group,donors=len(z),available=int(z[col].notna().sum()),
                missing_cases=','.join(map(str,z.loc[z[col].isna(),'case_id'])),minimum=z[col].min(),maximum=z[col].max()))
    pd.DataFrame(rows).to_csv(OUT/'glycemic_measure_availability.csv',index=False)
    methods=dict(original_h10=dict(status='Historical preapproval proposal',not_identified_in_supplied_paper_search=[
        'Censored predictive log-score comparison','Strongly regularized censored predictive regression',
        'Baseline-conditioned parametric benchmark']),
        proposed_typical_alternative=dict(status='Historical preapproval proposal; current specification is approved_analysis_specification.json',
            baseline='log C-peptide ~ sampled Beta density + age + admission glucose',
            augmented='baseline + one fixed equal-weight boundary-support / immune-access score',
            primary_population='19 donors with measured C-peptide; do not treat the censoring bound as observed',
            secondary_population='20 donors in separately labeled censor-bound substitution sensitivities',
            comparison='Nested ordinary regressions: added-score partial F test, adjusted R-squared; donor-held-out RMSE',
            glycemic_choice='Admission glucose covers all 20 donors and refers to a specified time; selected without outcome association checks',
            preprocessing='Training-only predictor/score scaling for held-out predictions; no feature selection',
            disease_sensitivity='Add disease indicators to both models; retain the same score definition'),
        local_precedents=[
            dict(path='/home/smith6jt/Documents/papers/Drotar2025.pdf',pdf_page=7,method='ANCOVA; outcome log transformation'),
            dict(path='/data/panc_CODEX/Islet-Explorer/literature/Heide.pdf',pdf_page=23,method='Linear regression fit and slope test'),
            dict(path='/home/smith6jt/Documents/papers/Williams2024.pdf',pdf_page=3,method='Random forest with separate 80/20 training/test sets; different model from proposed clinical regression'),
            dict(path='/home/smith6jt/Documents/papers/CELL-REPORTS-D-26-04511.pdf',pdf_page=12,method='RMSE for cell-count classifier validation; supplied manuscript draft, not the same clinical prediction task'),
            dict(path='/data/panc_CODEX/Islet-Explorer/literature/Richardson2023.pdf',pdf_page=3,method='Below-detection C-peptide assigned zero in that paper')],
        score_definition='0.5 * (training-z(mean_boundary_coverage) - training-z(mean_immune_localization)); equal weights and directions fixed before clinical fitting',
        measurement_definition='The score uses existing tissue measurements; the clinical outcome is recorded donor C-peptide')
    (OUT/'method_review.json').write_text(json.dumps(methods,indent=2)+'\n')
    plot=d[['case_id','donor_status','tissue_beta_density_mm2','C_peptide_ng_mL','C_peptide_censored','age_years_npod','admission_glucose_mg_dL','mean_boundary_coverage','mean_immune_localization']]
    plot.to_csv(OUT/'prism_clinical_input_points.csv',index=False)
    colors={'ND':'#4477AA','Aab+':'#CC6633','T1D':'#228833'}
    fig,ax=plt.subplots(figsize=(6,4.6))
    for group in ['ND','Aab+','T1D']:
        z=d[d.donor_status.eq(group)]
        measured=z[~z.C_peptide_censored]
        ax.scatter(measured.tissue_beta_density_mm2,measured.C_peptide_ng_mL,s=45,color=colors[group],label=group,
            edgecolor='white',linewidth=.5,zorder=3)
        for r in z[z.C_peptide_censored].itertuples():
            ax.scatter([r.tissue_beta_density_mm2],[r.C_peptide_ng_mL],s=60,color=colors[group],marker='v',zorder=3)
            ax.annotate(f'{r.case_id}: <0.02',(r.tissue_beta_density_mm2,r.C_peptide_ng_mL),xytext=(10,5),
                textcoords='offset points',fontsize=9,color=colors[group])
            ax.annotate('',xy=(r.tissue_beta_density_mm2,.012),xytext=(r.tissue_beta_density_mm2,.019),
                arrowprops=dict(arrowstyle='->',color=colors[group],lw=1))
    ax.set_yscale('log');ax.set_ylim(.01,65)
    ax.set_xlabel('Sampled β-cell quantity (cells/mm² tissue)')
    ax.set_ylabel('Recorded C-peptide (ng/mL)')
    ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(OUT/'clinical_quantity_overview.png',dpi=220);fig.savefig(OUT/'clinical_quantity_overview.svg')
    plt.close(fig)
    (OUT/'clinical_quantity_figure_legend.txt').write_text(
        'Recorded C-peptide and sampled Beta-cell quantity. Each point is one donor (ND n=7, Aab+ n=7, T1D n=6). '
        'Sampled quantity is the number of unique edge-eligible typed Beta cells inside the saved Tissue geometry divided by its area. '
        'C-peptide is displayed on a logarithmic axis. Nineteen circles represent measured concentrations; the downward triangle and arrow '
        'show donor 6550 recorded as below 0.02 ng/mL, rather than treating the bound as a measured value. '
        'This descriptive plot contains no model fit or hypothesis test. Tissue quantity is not whole-organ Beta-cell mass, and donor C-peptide was not a standardized stimulation assay.\n')
    (OUT/'data_measurement_review.txt').write_text('H10 measurement availability\n\n'
        'All 20 donors have sampled Beta quantity, age, admission glucose, and peak glucose. HbA1c is missing in donors 6090 and 6232. '
        'Nineteen donors have measured C-peptide; donor 6550 is recorded below 0.02 ng/mL.\n\n'
        'The method-search record is in method_review.json; the executed specification is approved_analysis_specification.json. '
        'The fixed score uses boundary coverage with a positive sign and immune core localization with a negative sign, with equal weights.\n')
    print('Saved clinical input plot, Prism table, glycemic availability, and exact method-review record; no clinical model fitted.')

if __name__=='__main__':main()
