#!/usr/bin/env python3
"""H9 donor-level coverage contrast, reusing the existing H2 coefficients."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h9-original-20261026-mpl')
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[2]
OUT=REPO/'output/h5_h10_original_plan_20261026/h9'
GROUPS=['ND','Aab+','T1D']
COLORS={'ND':'#4477AA','Aab+':'#CC6633','T1D':'#228833'}

def main():
    all_donors=pd.read_csv(OUT/'donor_inputs.csv')
    d=all_donors[all_donors.donor_status.isin(['Aab+','T1D'])].copy()
    objects=pd.read_csv(OUT/'analysis_objects_existing_h2.csv')
    assert len(d)==13 and 78 not in d.case_id.values
    beta_center=d.tissue_beta_density_mm2.mean()
    design=pd.DataFrame({'Intercept':np.ones(len(d)),'Aab_minus_T1D':d.donor_status.eq('Aab+').astype(float),
                         'beta_density_per_100':(d.tissue_beta_density_mm2-beta_center)/100},index=d.index)
    models={}
    model_rows=[]
    for name,extra in [('supporting_all_donors',[]),('age_adjusted_sensitivity',['age_years_npod']),
                       ('regional_immune_adjusted_sensitivity',['tissue_immune_density_mm2'])]:
        x=design.copy()
        for col in extra:
            x[col]=(d[col]-d[col].mean())/d[col].std()
        m=sm.OLS(d.coverage_core_immune_slope,x).fit()
        models[name]=m
        ci=m.conf_int().loc['Aab_minus_T1D']
        model_rows.append(dict(model=name,contrast='Aab+ minus T1D',estimate=m.params['Aab_minus_T1D'],
            SE=m.bse['Aab_minus_T1D'],ci_low=ci.iloc[0],ci_high=ci.iloc[1],t=m.tvalues['Aab_minus_T1D'],
            p_value=m.pvalues['Aab_minus_T1D'],residual_df=m.df_resid,donors=13,
            full_rank=np.linalg.matrix_rank(x)==x.shape[1],r_squared=m.rsquared,
            formula='coverage_core_immune_slope ~ Aab_vs_T1D + tissue_beta_density_mm2'+''.join(' + '+v for v in extra)))
        (OUT/f'{name}_model_summary.txt').write_text(m.summary().as_text()+'\n')
    results=pd.DataFrame(model_rows)
    results.to_csv(OUT/'adjusted_group_contrasts.csv',index=False)
    primary=models['supporting_all_donors']
    d['fitted_slope']=primary.fittedvalues
    d['residual']=primary.resid
    d['leverage']=primary.get_influence().hat_matrix_diag
    d.to_csv(OUT/'donor_model_values.csv',index=False)

    overlap=[]
    def support_rows(frame,variable,scope):
        a=frame.loc[frame.donor_status.eq('Aab+'),variable]
        t=frame.loc[frame.donor_status.eq('T1D'),variable]
        lower=max(a.min(),t.min()); upper=min(a.max(),t.max())
        overlap.append(dict(scope=scope,variable=variable,aab_min=a.min(),aab_max=a.max(),
            t1d_min=t.min(),t1d_max=t.max(),common_low=lower,common_high=upper,
            aab_values_in_common=int(a.between(lower,upper).sum()),t1d_values_in_common=int(t.between(lower,upper).sum())))
        return lower,upper
    low,high=support_rows(d,'tissue_beta_density_mm2','donors')
    common=d[d.tissue_beta_density_mm2.between(low,high)].copy()
    common.to_csv(OUT/'donors_in_common_beta_range.csv',index=False)
    for col in ['mean_boundary_coverage','mean_immune_localization','tissue_immune_density_mm2','age_years_npod']:
        support_rows(d,col,'donor_summaries_all')
    subset=objects[objects.case_id.isin(common.case_id)].copy()
    for col in ['coliv_coverage_2um','field_density','n_immune_nearby']:
        support_rows(subset,col,'objects_in_donors_with_common_beta_range')
    pd.DataFrame(overlap).to_csv(OUT/'overlap_ranges.csv',index=False)
    coverage_limits=(max(subset.loc[subset.donor_status.eq('Aab+'),'coliv_coverage_2um'].min(),subset.loc[subset.donor_status.eq('T1D'),'coliv_coverage_2um'].min()),
                     min(subset.loc[subset.donor_status.eq('Aab+'),'coliv_coverage_2um'].max(),subset.loc[subset.donor_status.eq('T1D'),'coliv_coverage_2um'].max()))
    field_limits=(max(subset.loc[subset.donor_status.eq('Aab+'),'field_density'].min(),subset.loc[subset.donor_status.eq('T1D'),'field_density'].min()),
                  min(subset.loc[subset.donor_status.eq('Aab+'),'field_density'].max(),subset.loc[subset.donor_status.eq('T1D'),'field_density'].max()))
    burden_limits=(max(subset.loc[subset.donor_status.eq('Aab+'),'n_immune_nearby'].min(),subset.loc[subset.donor_status.eq('T1D'),'n_immune_nearby'].min()),
                   min(subset.loc[subset.donor_status.eq('Aab+'),'n_immune_nearby'].max(),subset.loc[subset.donor_status.eq('T1D'),'n_immune_nearby'].max()))
    subset['within_joint_measured_ranges']=subset.coliv_coverage_2um.between(*coverage_limits)&subset.field_density.between(*field_limits)&subset.n_immune_nearby.between(*burden_limits)
    joint=subset.groupby(['case_id','donor_status']).agg(objects=('islet_id','size'),
        objects_in_joint_ranges=('within_joint_measured_ranges','sum'),core_immune_cells=('n_immune_core','sum')).reset_index()
    joint.to_csv(OUT/'joint_overlap_object_counts.csv',index=False)
    support=dict(beta_common_range=[low,high],donors_in_beta_common_range=common.case_id.tolist(),
        aab_donors_in_beta_common_range=int(common.donor_status.eq('Aab+').sum()),
        t1d_donors_in_beta_common_range=int(common.donor_status.eq('T1D').sum()),
        residual_df_for_three_parameter_model_in_common_range=int(len(common)-3),
        primary_original_support_status='Not estimable with residual uncertainty: insufficient donor replication in the shared Beta-quantity range',
        full_sample_model_status='Algebraically identifiable; depends on a shared linear Beta-quantity adjustment outside the narrow common range',
        source_coefficients='Existing H2 separate-donor adjusted slopes; exposure and continuous controls z within donor, immune-density outcome z globally',
        p_value_scope='All reported model p values are supporting extrapolation-dependent contrasts. The original overlap-supported primary test is unavailable and has no p value for the ten-hypothesis family.')
    (OUT/'support_assessment.json').write_text(json.dumps(support,indent=2)+'\n')

    fig,axs=plt.subplots(1,2,figsize=(10.4,4.5))
    ax=axs[0]
    arrays=[all_donors.loc[all_donors.donor_status.eq(g),'coverage_core_immune_slope'].values for g in GROUPS]
    boxes=ax.boxplot(arrays,positions=range(3),widths=.46,patch_artist=True,showfliers=False,
                     medianprops=dict(color='#222222',lw=1.4))
    for i,(g,values,b) in enumerate(zip(GROUPS,arrays,boxes['boxes'])):
        b.set_facecolor(COLORS[g]);b.set_alpha(.18);b.set_edgecolor(COLORS[g])
        ax.scatter(i+np.linspace(-.08,.08,len(values)),values,s=40,color=COLORS[g],edgecolor='white',lw=.5,zorder=3)
    ax.axhline(0,color='#888888',ls='--',lw=.8)
    ax.set_xticks(range(3),['ND\n(reference)','Aab+','T1D'])
    ax.set_xlabel('Disease state');ax.set_ylabel('Coverage–immune-density association\n(standardized adjusted coefficient)')
    ax=axs[1]
    ax.axvspan(low,high,color='#eeeeee',zorder=0)
    line_rows=[]
    for g in ['Aab+','T1D']:
        z=d[d.donor_status.eq(g)]
        ax.scatter(z.tissue_beta_density_mm2,z.coverage_core_immune_slope,s=43,color=COLORS[g],label=g,zorder=3)
        grid=np.linspace(z.tissue_beta_density_mm2.min(),z.tissue_beta_density_mm2.max(),100)
        xp=pd.DataFrame({'Intercept':1.,'Aab_minus_T1D':float(g=='Aab+'),'beta_density_per_100':(grid-beta_center)/100})
        yp=primary.predict(xp)
        ax.plot(grid,yp,color=COLORS[g],lw=1.4)
        line_rows.append(pd.DataFrame({'donor_status':g,'tissue_beta_density_mm2':grid,'fitted_slope':yp}))
        for r in z[z.case_id.isin(common.case_id)].itertuples():
            ax.annotate(str(r.case_id),(r.tissue_beta_density_mm2,r.coverage_core_immune_slope),
                xytext=(4,5),textcoords='offset points',fontsize=8,color=COLORS[g])
    ax.axhline(0,color='#888888',ls='--',lw=.8)
    ax.set_xlabel('Sampled β-cell quantity (cells/mm² tissue)')
    ax.set_ylabel('Coverage–immune-density association\n(standardized adjusted coefficient)')
    ax.legend(frameon=False,fontsize=9)
    for a in axs:a.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT/'coverage_association_and_beta_quantity.png',dpi=220)
    fig.savefig(OUT/'coverage_association_and_beta_quantity.svg')
    plt.close(fig)
    all_donors[['case_id','donor_status','coverage_core_immune_slope','tissue_beta_density_mm2']].to_csv(OUT/'prism_donor_points.csv',index=False)
    pd.concat(line_rows).to_csv(OUT/'prism_adjusted_regression_lines.csv',index=False)
    r=results.iloc[0]
    legend=("Boundary coverage, immune localization and sampled Beta-cell quantity. "
        "Each point is one donor (ND, n=7; Aab+, n=7; T1D, n=6). The left panel shows the existing adjusted "
        "within-donor regression coefficient linking nonvascular ColIV-positive coverage in a 2-µm boundary band "
        "to core immune density. Exposure and continuous controls were standardized within donor; density was "
        "standardized across all eligible objects. Controls were peri immune density, core area, perimeter, tissue-edge "
        "distance, capillary area fraction, and object class. Negative coefficients indicate lower core immune density "
        "with greater boundary coverage. Boxes show median and interquartile range; whiskers extend to observations "
        "within 1.5 times that range; all donor points are shown. ND is a descriptive reference. The right panel relates "
        "the Aab+/T1D donor coefficients to unique edge-eligible typed Beta cells per mm² of sampled Tissue. Lines "
        "show the supporting all-donor linear model with disease and Beta quantity (13 donors; 10 residual degrees of freedom). "
        f"The adjusted Aab+-minus-T1D contrast is {r.estimate:+.4f} (95% CI {r.ci_low:+.4f} to {r.ci_high:+.4f}; "
        f"two-sided t test, p={r.p_value:.3f}). The shaded shared Beta-quantity range ({low:.2f}–{high:.2f} cells/mm²) "
        "contains only one Aab+ and two T1D donors, identified by case number. The adjusted comparison therefore "
        "depends on a linear adjustment beyond this narrow overlap; a three-parameter model restricted to the shared "
        "range has no residual degrees of freedom. Therefore the original plan's overlap-supported primary contrast is not estimable with residual uncertainty. "
        "The displayed all-donor model is a separate supplementary calculation. No significance symbol is shown because p is greater than 0.05. No biological branch judgment is assigned.")
    (OUT/'figure_legend.txt').write_text(legend+'\n')
    report=["H9: boundary-coverage association adjusted for donor-wide Beta quantity",'',
        "Original-plan primary status: NOT ESTIMABLE with residual uncertainty within the shared Beta-quantity range.",'',
        f"Supporting full 13-donor Aab+-minus-T1D contrast: {r.estimate:+.6f}, 95% CI [{r.ci_low:+.6f}, {r.ci_high:+.6f}], p={r.p_value:.6f}, residual df={int(r.residual_df)}.",
        "The original hypothesis specified a negative Aab+-minus-T1D contrast. No biological interpretation is assigned here.",
        '',f"The shared Beta-quantity range is {low:.6f} to {high:.6f} cells/mm²: one Aab+ donor and two T1D donors. "
        "All thirteen donors have core immune events, and object-level coverage/exposure ranges overlap, but these objects do not increase the number of independent donors in shared Beta support. "
        "The original overlap-supported conditional contrast is not estimable with residual uncertainty from those three donors. The supplementary all-donor adjusted model includes donors outside this shared range.",
        '',"Sensitivity contrasts:"]
    for row in results.iloc[1:].itertuples():report.append(f"{row.model}: {row.estimate:+.6f} [{row.ci_low:+.6f}, {row.ci_high:+.6f}], p={row.p_value:.6f}.")
    report+=['',"Existing corrected counts, H2 threshold measurements, and H2 within-donor coefficients were reused. The added Beta census counts each valid cell once within saved Tissue geometry; it uses neither duplicated peri memberships nor only immune-exposed objects.",
        '',"Methods in supplied papers: Drotar2025.pdf, PDF p7 (ANCOVA with a donor covariate); Barlow2025.pdf, PDF p22 (within-donor predictor standardization); Bodenmiller2025.pdf, pp17–18 (donor dot/boxplots). The specific donor-slope endpoint follows the original H9 plan.",
        '',"The audit records the original plan's overlap requirement and model estimability. Result interpretation and biological branch assignment are left unassigned under the owner's instruction."]
    (OUT/'RESULTS.txt').write_text('\n'.join(report)+'\n')
    branch=dict(hypothesis='H9',original_primary_status='not_estimable',primary_p_value=None,
        if_supported_activated=None,if_meaningful_difference_excluded_activated=None,
        reason='Only one Aab+ and two T1D donors are in the shared Beta-quantity range; the three-parameter adjusted model has zero residual degrees of freedom.',
        supporting_all_donor_p=float(r.p_value),
        branch_assignment='Not assigned: owner requested numerical reporting and plan-compliance checking without autonomous interpretation.')
    (OUT/'conditional_branch.json').write_text(json.dumps(branch,indent=2)+'\n')
    print('\n'.join(report))

if __name__=='__main__':main()
