#!/usr/bin/env python3
"""Approved donor H9/H10 models and existing-app profiles; no INS endpoints."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/app_h9_h10_mpl')
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h10_approved_models_20261026 as clinical

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'output/app_hypotheses_20261026/h9_h10'
OLD = ROOT/'output/h5_h10_original_plan_20261026'
APP = ROOT/'data/development/uniform-20donor-20261022-v10'
SHARED = ROOT/'output/app_hypotheses_20261026/shared/object_measurements.parquet'
H2_OBJECTS = ROOT/'output/app_hypotheses_20261026/h1_h4/data/h2_analysis_objects.csv'
H2_SLOPES = ROOT/'output/app_hypotheses_20261026/h1_h4/statistics/h2_donor_values.csv'
GROUPS = ['ND', 'Aab+', 'T1D']
COLORS = {'ND':'#4477AA', 'Aab+':'#CC6633', 'T1D':'#228833'}
MANIFEST = []


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()


def save(fig, folder, stem, legend, sources):
    fig.tight_layout()
    for ext in ['png','pdf','svg']:
        fig.savefig(folder/'figures'/f'{stem}.{ext}', dpi=220, bbox_inches='tight')
    plt.close(fig)
    (folder/'figures'/f'{stem}_legend.txt').write_text(legend+'\n')
    MANIFEST.append(dict(family=folder.name, stem=stem,
        figure_files=[str((folder/'figures'/f'{stem}.{ext}').relative_to(OUT)) for ext in ['png','pdf','svg']],
        legend=str((folder/'figures'/f'{stem}_legend.txt').relative_to(OUT)),
        source_csv=[str((folder/'data'/s).relative_to(OUT)) for s in sources]))


def prepare():
    for f in ['h9','h10','profiles']:
        for s in ['data','figures']:(OUT/f/s).mkdir(parents=True,exist_ok=True)
    # Existing tissue census is rebound to its exact retained cell/geometry inputs.
    prior=json.loads((OLD/'h9/source_provenance.json').read_text())
    audit=[]
    for path, expected in prior['source_files'].items():
        current=sha(path)
        audit.append(dict(path=path, previous_sha256=expected, current_sha256=current, unchanged=current==expected))
    assert all(a['unchanged'] for a in audit), 'A source changed; refresh the affected measurement before reuse.'
    pd.DataFrame(audit).to_csv(OUT/'source_binding_audit.csv',index=False)
    d=pd.read_csv(OLD/'h10/donor_inputs.csv').sort_values('case_id').reset_index(drop=True)
    meta=pd.read_csv(APP/'donor_metadata.csv').sort_values('case_id').reset_index(drop=True)
    for c in meta.columns:
        pd.testing.assert_series_equal(d[c],meta[c],check_dtype=False,check_names=False)
    objects=pd.read_csv(H2_OBJECTS,low_memory=False)
    slopes=pd.read_csv(H2_SLOPES)
    check=objects.groupby('case_id').agg(mean_boundary_coverage=('coliv_coverage_2um','mean'),
        mean_immune_localization=('immune_localization_fraction','mean')).reindex(d.case_id)
    np.testing.assert_allclose(check.to_numpy(),d[check.columns].to_numpy(),rtol=1e-12,atol=1e-14)
    np.testing.assert_allclose(slopes.set_index('case_id').loc[d.case_id,'standardized_beta'],d.coverage_core_immune_slope)
    assert len(d)==20 and d.case_id.nunique()==20 and 78 not in d.case_id.values
    assert d.groupby('donor_status').size().to_dict()=={'ND':7,'Aab+':7,'T1D':6}
    assert d[clinical.BASE+clinical.COMPONENTS+['C_peptide_ng_mL']].notna().all().all()
    assert d.loc[d.C_peptide_censored,'case_id'].tolist()==[6550]
    for family in ['h9','h10']:
        d.to_csv(OUT/family/'data/donor_inputs.csv',index=False)
        pd.read_csv(OLD/'h9/tissue_beta_quantity.csv').to_csv(OUT/family/'data/tissue_beta_quantity.csv',index=False)
    objects.to_csv(OUT/'h9/data/analysis_objects.csv',index=False)
    slopes.to_csv(OUT/'h9/data/source_donor_coefficients.csv',index=False)
    provenance=dict(sources={str(p.relative_to(ROOT)):sha(p) for p in [SHARED,H2_OBJECTS,H2_SLOPES,APP/'donor_metadata.csv',Path(__file__)]},
        original_tissue_census=prior['quantity_definition'],source_bindings_unchanged=True,
        structural_area_correction='The fixed score contains no summed annotation area and no recovered structural-cell counts. H2 already uses hf.union_local for capillary masks before nonvascular boundary coverage and capillary core-area fraction.',
        evidence='scripts/research/h2_coliv_boundary_20261025.py:74-85',
        score_components_recomputed_from_current_bound_h2_objects_and_identical=True,
        structural_sidecars='Included in shared provenance; app phenotype labels are not overwritten.',
        original_quantitative_INS_followup='Omitted by current user instruction; no replacement tissue outcome is fitted.')
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    return d,objects


def h9(d,objects):
    folder=OUT/'h9';data=folder/'data'
    z=d[d.donor_status.isin(['Aab+','T1D'])].copy()
    bmean=z.tissue_beta_density_mm2.mean()
    def design(q, extra=None):
        x=pd.DataFrame({'Intercept':1., 'Aab_minus_T1D':q.donor_status.eq('Aab+').astype(float),
            'beta_density_per_100':(q.tissue_beta_density_mm2-bmean)/100},index=q.index)
        if extra:x[extra]=(q[extra]-z[extra].mean())/z[extra].std()
        return x
    models={};rows=[]
    for label,extra in [('all_13_donors',None),('age_adjusted','age_years_npod'),('regional_immune_adjusted','tissue_immune_density_mm2')]:
        x=design(z,extra);m=sm.OLS(z.coverage_core_immune_slope,x).fit();models[label]=m
        ci=m.conf_int().loc['Aab_minus_T1D']
        rows.append(dict(analysis=label,comparison='Aab+ minus T1D',estimate=m.params.Aab_minus_T1D,
            lower_95=ci.iloc[0],upper_95=ci.iloc[1],SE=m.bse.Aab_minus_T1D,p_value=m.pvalues.Aab_minus_T1D,
            statistic=m.tvalues.Aab_minus_T1D,residual_df=m.df_resid,n_donors=len(z),
            scope='Supplementary; includes donors outside shared Beta-quantity range',
            formula='H2 within-donor coefficient ~ disease + tissue Beta density'+(' + '+extra if extra else '')))
        (data/f'{label}_model_summary.txt').write_text(m.summary().as_text()+'\n')
    results=pd.DataFrame(rows);results.to_csv(data/'adjusted_group_contrasts.csv',index=False)
    ranges=[]
    def overlap(q,col,scope):
        a=q.loc[q.donor_status.eq('Aab+'),col];t=q.loc[q.donor_status.eq('T1D'),col]
        lo=max(a.min(),t.min());hi=min(a.max(),t.max())
        ranges.append(dict(scope=scope,variable=col,aab_min=a.min(),aab_max=a.max(),t1d_min=t.min(),t1d_max=t.max(),
            common_low=lo,common_high=hi,aab_values_in_common=int(a.between(lo,hi).sum()),t1d_values_in_common=int(t.between(lo,hi).sum())))
        return lo,hi
    lo,hi=overlap(z,'tissue_beta_density_mm2','donors')
    z['within_shared_beta_range']=z.tissue_beta_density_mm2.between(lo,hi)
    common=z[z.within_shared_beta_range].copy();common.to_csv(data/'donors_in_common_beta_range.csv',index=False)
    for c in ['mean_boundary_coverage','mean_immune_localization','tissue_immune_density_mm2','age_years_npod']:overlap(z,c,'donors')
    subset=objects[objects.case_id.isin(common.case_id)].copy();within=np.ones(len(subset),bool)
    for c in ['coliv_coverage_2um','field_density','n_immune_nearby']:
        bounds=overlap(subset,c,'objects_of_donors_in_shared_beta_range');within &= subset[c].between(*bounds)
    subset['within_joint_measured_ranges']=within
    subset.to_csv(data/'objects_in_common_beta_range.csv',index=False)
    counts=subset.groupby(['case_id','donor_status']).agg(n_objects=('islet_id','size'),n_in_joint_ranges=('within_joint_measured_ranges','sum'),
        core_immune_cells=('n_immune_core','sum')).reset_index();counts.to_csv(data/'joint_overlap_counts.csv',index=False)
    pd.DataFrame(ranges).to_csv(data/'overlap_ranges.csv',index=False)
    residual_df=len(common)-np.linalg.matrix_rank(design(common))
    primary=dict(hypothesis='H9',comparison='Aab+ versus T1D adjusted H2 slope within shared Beta-quantity range',
        test='Donor-level linear-model group coefficient, adjusted for tissue Beta density',p_value=np.nan,
        n_donors=len(common),n_objects=len(subset),status='unavailable_zero_residual_df')
    if residual_df>0 and np.linalg.matrix_rank(design(common))==3:
        fit=sm.OLS(common.coverage_core_immune_slope,design(common)).fit()
        primary.update(p_value=fit.pvalues.Aab_minus_T1D,status='computed')
        (data/'shared_range_model_summary.txt').write_text(fit.summary().as_text()+'\n')
    support=dict(beta_quantity_shared_range=[lo,hi],shared_range_case_ids=common.case_id.tolist(),
        shared_range_group_n=common.groupby('donor_status').size().to_dict(),residual_df=int(residual_df),
        n_objects_in_shared_donors=len(subset),n_objects_in_joint_ranges=int(within.sum()),
        primary_status=primary['status'],branch_decision='unassigned_user_judgment')
    (folder/'technical_estimability.json').write_text(json.dumps(support,indent=2)+'\n')
    d[['case_id','donor_status','coverage_core_immune_slope','tissue_beta_density_mm2']].to_csv(data/'donor_slopes.csv',index=False)
    z.to_csv(data/'conditional_donor_values.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    for j,g in enumerate(GROUPS):
        v=d.loc[d.donor_status.eq(g),'coverage_core_immune_slope'].to_numpy()
        axes[0].boxplot(v,positions=[j],widths=.43,showfliers=False,medianprops={'color':'black'})
        axes[0].scatter(j+np.linspace(-.08,.08,len(v)),v,s=37,color=COLORS[g],zorder=3)
    axes[0].set(xticks=range(3),xticklabels=GROUPS,xlabel='Disease state',ylabel='Adjusted coverage–immune-density coefficient')
    axes[1].axvspan(lo,hi,color='#dddddd',alpha=.5)
    lines=[]
    for g in ['Aab+','T1D']:
        q=z[z.donor_status.eq(g)]
        axes[1].scatter(q.tissue_beta_density_mm2,q.coverage_core_immune_slope,s=40,color=COLORS[g],label=g)
        grid=pd.DataFrame({'tissue_beta_density_mm2':np.linspace(q.tissue_beta_density_mm2.min(),q.tissue_beta_density_mm2.max(),100),'donor_status':g})
        grid['fitted_coefficient']=models['all_13_donors'].predict(design(grid));lines.append(grid)
        axes[1].plot(grid.tissue_beta_density_mm2,grid.fitted_coefficient,color=COLORS[g])
        for r in q[q.within_shared_beta_range].itertuples():axes[1].annotate(str(r.case_id),(r.tissue_beta_density_mm2,r.coverage_core_immune_slope),xytext=(4,4),textcoords='offset points',fontsize=8)
    pd.concat(lines).to_csv(data/'conditional_regression_lines.csv',index=False)
    axes[1].set(xlabel='Sampled β-cell density (cells/mm² tissue)',ylabel='Adjusted coverage–immune-density coefficient');axes[1].legend(frameon=False)
    for ax in axes:ax.axhline(0,lw=.8,ls='--',color='#888888')
    r=results.iloc[0]
    legend=('Each point is one donor (ND 7, Aab+ 7, T1D 6); ND is shown as reference. The coefficient is the existing H2 adjusted within-donor association of 2-µm nonvascular ColIV boundary coverage with core immune density. Exposure and continuous controls are standardized within donor; the immune-density response is standardized globally. Controls: peri immune density, core area, perimeter, tissue-edge distance, capillary core-area fraction and object class. Boxes show median/IQR, whiskers extend to 1.5 IQR, and all donor points are displayed. The right panel contains 13 Aab+/T1D donors; lines are the supplementary linear model of the donor coefficient on disease and unique tissue Beta-cell density. '
        f'Aab+ minus T1D = {r.estimate:.6f}, 95% CI [{r.lower_95:.6f}, {r.upper_95:.6f}], two-sided t p={r.p_value:.6g}, residual df={int(r.residual_df)}. The shaded shared Beta-density range is [{lo:.6f}, {hi:.6f}] cells/mm²; labeled donors are inside that range. The restricted three-parameter design has {residual_df} residual degrees of freedom. No p value is assigned to the original shared-range contrast when its residual degrees of freedom are zero. No significance asterisk is shown for an unavailable primary test; the supplementary contrast is reported separately.')
    save(fig,folder,'coverage_association_and_beta_quantity',legend,['donor_slopes.csv','conditional_donor_values.csv','conditional_regression_lines.csv','adjusted_group_contrasts.csv'])
    return primary,support,results


def h10(d):
    folder=OUT/'h10';data=folder/'data';clinical.OUT=data
    d=d.copy();d['log_c_peptide']=np.log(d.C_peptide_ng_mL)
    measured=d[~d.C_peptide_censored].copy()
    comparisons=[]
    p,pred,fits=clinical.ordinary_analysis(measured,'log_c_peptide','clinical',False);comparisons.append(p)
    c,_,_=clinical.ordinary_analysis(measured,'log_c_peptide','clinical',True);comparisons.append(c)
    pd.DataFrame(comparisons).to_csv(data/'ordinary_model_comparisons.csv',index=False)
    cens=d[d.C_peptide_censored].copy();cp=cens[['case_id','donor_status','C_peptide_ng_mL','C_peptide_censored']].copy()
    for aug in [False,True]:
        name='augmented' if aug else 'baseline';x,_=clinical.design(measured,cens,aug,False)
        cp[name+'_prediction']=fits[name].predict(x).to_numpy()
    cp['prediction_basis']='Trained on 19 measured donors; censored outcome unscored'
    cp.to_csv(data/'censored_donor_unscored_primary_predictions.csv',index=False)
    t=d[d.donor_status.eq('T1D')][['case_id','donor_status','C_peptide_ng_mL','C_peptide_censored']].merge(
        pd.concat([pred[['case_id','baseline_prediction','augmented_prediction']],cp[['case_id','baseline_prediction','augmented_prediction']]]),on='case_id',validate='one_to_one')
    for name in ['baseline','augmented']:
        t[name+'_predicted_median_ng_mL']=np.exp(t[name+'_prediction'])
        t[name+'_log_error']=np.where(t.C_peptide_censored,np.nan,t[name+'_prediction']-np.log(t.C_peptide_ng_mL))
    t['prediction_basis']=np.where(t.C_peptide_censored,'Trained on 19 measured donors; outcome unscored','Held out; trained on 18 measured donors')
    t.to_csv(data/'all_six_T1D_observations_and_predictions.csv',index=False)
    coefs=[];preds=[];diags=[];lrtests=[]
    for disease in [False,True]:
        ds=[]
        for aug in [False,True]:
            rows,vals,diag=clinical.censored_fit(d,aug,disease);coefs+=rows;preds.append(vals);diags.append(diag);ds.append(diag)
        lr=2*(ds[1]['log_likelihood']-ds[0]['log_likelihood']);assert lr>=-1e-8
        lrtests.append(dict(analysis=ds[0]['analysis'],n_donors=20,censored_donors=1,test='Likelihood-ratio test',statistic=max(0.,lr),df=1,p_value=stats.chi2.sf(max(0.,lr),1)))
    pd.DataFrame(coefs).to_csv(data/'censored_model_coefficients.csv',index=False)
    pd.concat(preds).to_csv(data/'censored_model_fitted_values.csv',index=False)
    pd.DataFrame(diags).to_csv(data/'censored_model_diagnostics.csv',index=False)
    pd.DataFrame(lrtests).to_csv(data/'censored_likelihood_comparisons.csv',index=False)
    x,scale=clinical.design(d,d,True,False)
    display=d[['case_id','donor_status','mean_boundary_coverage','mean_immune_localization','tissue_beta_density_mm2','C_peptide_ng_mL','C_peptide_censored']].copy()
    display['spatial_score_all20_display']=x.spatial_score
    display.to_csv(data/'fixed_score_and_clinical_measurements.csv',index=False)
    scale.to_csv(data/'fixed_score_all20_display_scaling.csv',index=False)
    specification=dict(status='Approved method reused with audited unchanged inputs',baseline=clinical.BASE,
        score='0.5 * (z(mean_boundary_coverage) - z(mean_immune_localization))',
        score_aggregation='Arithmetic donor mean over existing eligible immune-exposed H2 objects; fixed components and equal weights',
        scaling='Full model uses its included donors; each held-out fold uses only training donors; all-20 score is for donor profile display and censored sensitivity',
        primary='19 measured donors; nested ordinary regression of natural log C-peptide; partial F test for added score',
        prediction='Leave one donor out; RMSE in natural-log C-peptide units, separate from the F test',
        censored_sensitivity='20 donors; normal log-density for measured values; normal log-CDF at log(0.02) for case6550; LR test',
        disease_sensitivity='Aab+ and T1D indicators in both models with ND reference',
        omitted_followup='Quantitative within-Beta INS follow-up removed under current user instruction; no replacement outcome fitted',
        cluster_features_in_clinical_models=False,imputation=False,feature_selection=False,branch_decision='unassigned_user_judgment')
    (folder/'analysis_specification.json').write_text(json.dumps(specification,indent=2)+'\n')
    return dict(hypothesis='H10',comparison='Addition of fixed spatial score beyond Beta quantity, age and admission glucose',
        test='Nested ordinary-regression partial F test',p_value=p['p_value'],n_donors=19,n_objects=np.nan,status='computed'),pd.DataFrame(comparisons),pd.DataFrame(lrtests)


def clinical_plots(comparisons,lrtests):
    folder=OUT/'h10';data=folder/'data'
    definitions=('Baseline predictors are unique sampled tissue Beta-cell density, age and admission glucose. The augmented model adds the fixed score 0.5 × [standardized mean nonvascular ColIV boundary coverage − standardized mean immune core localization]. Its components and equal weights are unchanged. ')
    common=('Each point is one donor. '+definitions+'Scaling is fitted on training donors only in every held-out fold, with no imputation, feature selection or tuning. ')
    for label,title in [('clinical',''),('clinical_disease_conditioned','Both models also include disease indicators (ND reference). ')]:
        q=pd.read_csv(data/f'{label}_held_out_predictions.csv');r=comparisons.set_index('analysis').loc[label]
        fig,axes=plt.subplots(1,2,figsize=(9.5,4.4),sharex=True,sharey=True)
        v=q[['observed','baseline_prediction','augmented_prediction']].to_numpy();lo=np.min(v);hi=np.max(v);pad=.08*(hi-lo)
        for ax,name,heading in zip(axes,['baseline','augmented'],['Baseline','With spatial score']):
            for g in GROUPS:
                a=q[q.donor_status.eq(g)];ax.scatter(a.observed,a[name+'_prediction'],color=COLORS[g],s=42,label=g)
            ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],ls='--',color='#888888',lw=.8)
            ax.set(xlabel='Observed ln C-peptide',title=heading,xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad))
        axes[0].set_ylabel('Held-out predicted ln C-peptide');axes[0].legend(frameon=False)
        legend=common+title+f'Nineteen measured donors (ND 7, Aab+ 7, T1D 5); outcome is ln(C-peptide in ng/mL). The dashed line is identity. Baseline/augmented held-out RMSE: {r.baseline_held_out_rmse:.6f}/{r.augmented_held_out_rmse:.6f} natural-log units. Full-sample added-score F(1,{int(r.df_denominator)})={r.statistic:.6f}, p={r.p_value:.6g}; adjusted R² {r.baseline_adjusted_r_squared:.6f}/{r.augmented_adjusted_r_squared:.6f}. F tests and held-out error are separate statistics; no test of RMSE change is performed.'
        save(fig,folder,label+'_held_out_prediction',legend,[label+'_held_out_predictions.csv','ordinary_model_comparisons.csv'])
    t=pd.read_csv(data/'all_six_T1D_observations_and_predictions.csv').sort_values('case_id');fig,ax=plt.subplots(figsize=(9,4.5));ix=np.arange(6);m=~t.C_peptide_censored
    ax.scatter(ix[m]-.17,t.loc[m,'C_peptide_ng_mL'],c='#222222',s=45,label='Measured')
    for off,name,col,marker in [(0,'baseline','#4477AA','s'),(.17,'augmented','#CC6633','D')]:ax.scatter(ix+off,t[name+'_predicted_median_ng_mL'],c=col,s=40,marker=marker,label='Baseline prediction' if name=='baseline' else 'With spatial score')
    for i,r in enumerate(t.itertuples()):
        if r.C_peptide_censored:
            ax.scatter(i-.17,.02,marker='v',facecolors='none',edgecolors='#222222',s=65,label='Below 0.02')
            ax.annotate('',xy=(i-.17,.011),xytext=(i-.17,.019),arrowprops=dict(arrowstyle='->',color='#222222'))
    ax.set(xticks=ix,xticklabels=t.case_id.astype(str),xlabel='T1D donor',ylabel='C-peptide (ng/mL)',yscale='log',ylim=(.008,7));ax.legend(frameon=False,ncol=2,fontsize=8)
    save(fig,folder,'individual_T1D_predictions','Each symbol is a recorded measurement, censoring bound or model prediction for the donor on the horizontal axis. '+definitions+'All six T1D donors are displayed. Five measured donors have predictions from the other 18 measured donors; all predictor and score-component scaling uses those training donors. The open downward triangle is the recorded <0.02 bound for donor 6550; it is not a measured concentration. Its model and scaling use the 19 measured donors and no squared error is assigned. Exponentiated predicted log values are predicted medians in ng/mL. The vertical axis is logarithmic. No test is displayed.', ['all_six_T1D_observations_and_predictions.csv'])
    q=pd.read_csv(data/'fixed_score_and_clinical_measurements.csv');fig,axes=plt.subplots(1,2,figsize=(10,4.4))
    for i,g in enumerate(GROUPS):
        a=q[q.donor_status.eq(g)]
        axes[0].scatter(i+np.linspace(-.08,.08,len(a)),a.spatial_score_all20_display,color=COLORS[g],s=40)
        for r in a.itertuples():
            axes[1].scatter(r.spatial_score_all20_display,r.C_peptide_ng_mL,marker='v' if r.C_peptide_censored else 'o',facecolors='none' if r.C_peptide_censored else COLORS[g],edgecolors=COLORS[g],s=45)
            if r.C_peptide_censored:axes[1].annotate('6550: <0.02',(r.spatial_score_all20_display,.02),xytext=(5,4),textcoords='offset points',fontsize=8)
    axes[0].set(xticks=range(3),xticklabels=GROUPS,xlabel='Disease state',ylabel='Fixed spatial score');axes[1].set(xlabel='Fixed spatial score',ylabel='Recorded C-peptide (ng/mL)',yscale='log')
    save(fig,folder,'fixed_spatial_score_and_c_peptide','Each point is one donor. '+definitions+'All 20 donors (7 ND, 7 Aab+, 6 T1D) are plotted. For this display, both score components are standardized across all 20 donors. This display scaling does not enter held-out fits. The open triangle marks donor 6550 at its recorded upper bound of 0.02 ng/mL; no measured value is assigned. The concentration axis is logarithmic. No unadjusted correlation or group test is performed on these display points.', ['fixed_score_and_clinical_measurements.csv','fixed_score_all20_display_scaling.csv'])
    co=pd.read_csv(data/'censored_model_coefficients.csv');q=co[co.model.eq('augmented')&co.term.eq('spatial_score')].copy();q.to_csv(data/'plot_censored_score_coefficients.csv',index=False)
    fig,ax=plt.subplots(figsize=(6.5,4.4));ix=np.arange(len(q));ax.errorbar(ix,q.estimate,yerr=np.array([q.estimate-q.lower_95,q.upper_95-q.estimate]),fmt='o',capsize=5,color='#333333')
    ax.axhline(0,color='#888888',ls='--',lw=.8);ax.set(xticks=ix,xticklabels=['Clinical covariates','Also adjusted for\ndisease state'],ylabel='Spatial-score coefficient for ln C-peptide\nwith 95% CI',xlim=(-.5,1.5))
    tests='; '.join(f"{r.analysis}: LR χ²(1)={r.statistic:.6f}, p={r.p_value:.6g}" for r in lrtests.itertuples())
    save(fig,folder,'censored_score_coefficients','Each point is the fitted spatial-score regression coefficient from one model; error bars are observed-information Wald 95% confidence intervals. '+definitions+'Both models use all 20 donors to fit the regression and standardize the continuous predictors and score components. The second model also includes disease indicators with ND as reference. Measured C-peptide values contribute normal log-densities; donor 6550 contributes the normal log-CDF below ln(0.02). No bound substitution is used. '+tests+'. These are full-sample sensitivity models, not held-out prediction scores. No significance symbols are shown because the sensitivity p values are not below 0.05.', ['plot_censored_score_coefficients.csv','censored_likelihood_comparisons.csv'])


def profiles(d,objects):
    folder=OUT/'profiles';data=folder/'data';allq=pd.read_parquet(SHARED);q=allq[allq.eligible_plotted].copy();q.case_id=q.case_id.astype(int)
    assert len(q)==1126 and q.case_id.nunique()==20
    ordered=d.assign(group=d.donor_status.map({g:i for i,g in enumerate(GROUPS)})).sort_values(['group','case_id'])
    denom=q.groupby('case_id').size();occ=[]
    for r in ordered.itertuples():
        for cl in range(10):
            n=int(((q.case_id==r.case_id)&(q['leiden_0.5']==cl)).sum())
            occ.append(dict(case_id=r.case_id,donor_status=r.donor_status,cluster=cl,n_objects=n,denominator_objects=int(denom[r.case_id]),occupancy_percent=100*n/denom[r.case_id]))
    occ=pd.DataFrame(occ);occ.to_csv(data/'donor_cluster_occupancy.csv',index=False)
    assert np.allclose(occ.groupby('case_id').occupancy_percent.sum(),100)
    reference=pd.read_csv(ROOT/'output/app_hypotheses_20261026/h5/data/donor_cluster_occupancy.csv')
    keys=['case_id','donor_status','cluster']
    chk=occ.merge(reference,on=keys,validate='one_to_one',suffixes=('_here','_shared'))
    assert len(chk)==200 and np.allclose(chk.occupancy_percent_here,chk.occupancy_percent_shared)
    fig,ax=plt.subplots(figsize=(8.6,8));mat=occ.pivot(index='case_id',columns='cluster',values='occupancy_percent').loc[ordered.case_id]
    im=ax.imshow(mat,cmap='viridis',vmin=0,vmax=100,aspect='auto');ax.set(xticks=range(10),xlabel='Saved composition cluster (Leiden 0.5)',yticks=range(20),yticklabels=[f'{r.donor_status}  {r.case_id}' for r in ordered.itertuples()]);fig.colorbar(im,ax=ax,label='Donor objects in cluster (%)')
    for line in [6.5,13.5]:ax.axhline(line,color='white',lw=1.5)
    pop=q.groupby('donor_status').size().to_dict()
    intro=f'The matched display cohort contains1126 objects (ND {pop["ND"]}, Aab+ {pop["Aab+"]}, T1D {pop["T1D"]}) from20 donors (7/7/6), with complete core and20µm peri regions and an existing composition embedding. '
    save(fig,folder,'donor_composition_cluster_occupancy',intro+'Each row is one donor. Percent denominator is that donor’s eligible embedded objects; all10 saved Leiden0.5 clusters are included, with zero for absent clusters. The embedding and labels are reused without refitting. Cluster occupancy is a donor profile and does not enter either clinical model. No hypothesis test is performed.', ['donor_cluster_occupancy.csv'])
    q['pt_bin']=np.minimum((q.pt_rank*20).astype(int),19);q['pt_midpoint']=(q.pt_bin+.5)/20
    features=['immune_density_core_mm2','immune_density_peri_mm2','immune_localization_fraction']
    keep=['islet_id','case_id','donor_status','leiden_0.5','pt_rank','pt_bin','pt_midpoint',*features,'n_immune_core','n_immune_peri','n_immune_nearby','core_area_um2','peri_area_um2']
    q[keep].to_csv(data/'core_peri_profile_objects.csv',index=False)
    bins=q.groupby(['case_id','donor_status','pt_bin','pt_midpoint'])[features].agg(['mean','count']).reset_index()
    bins.columns=['_'.join(c).rstrip('_') if isinstance(c,tuple) else c for c in bins.columns];bins.to_csv(data/'donor_pseudotime_core_peri_profiles.csv',index=False)
    fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True,sharey='row')
    labels=['Core immune density (cells/mm²)','Peri immune density (cells/mm²)','Core immune fraction of nearby immune cells']
    for j,(f,label) in enumerate(zip(features,labels)):
        for k,g in enumerate(GROUPS):
            for case,a in bins[bins.donor_status.eq(g)].groupby('case_id'):
                axes[j,k].plot(a.pt_midpoint,a[f+'_mean'],'o-',lw=.8,ms=2,alpha=.65,color=COLORS[g])
            axes[j,k].set_title(g if j==0 else '');axes[j,k].set_xlim(0,1)
            if k==0:axes[j,k].set_ylabel(label)
            if j==2:axes[j,k].set_xlabel('Pseudotime (rank-normalized)')
    save(fig,folder,'donor_pseudotime_core_peri_profiles',intro+'Each line is one donor; each point is the arithmetic mean of observed object measurements in one of20 fixed app pseudotime rank bins. Bins without observations are not assigned zero. Core density=count/core area; peri density=count/assessed peri area. The20µm peri band is tissue-clipped and excludes Islet/Cluster cores. A cell may be associated with more than one peri object; these are object-associated profiles, not a unique tissue census. Localization is core/(core+peri) immune count, undefined if no nearby immune cells. Pseudotime values and cohort ranks are reused without refitting or reranking. No hypothesis tests are performed.', ['core_peri_profile_objects.csv','donor_pseudotime_core_peri_profiles.csv'])
    # Keep all original H2 eligible donors, independently of representation eligibility.
    full=objects.copy();full['immune_density_core_mm2']=full.n_immune_core*1e6/full.core_area_um2;full['immune_density_peri_mm2']=full.n_immune_peri*1e6/full.peri_area_um2
    agg=full.groupby(['case_id','donor_status']).agg(n_objects=('islet_id','size'),core_density=('immune_density_core_mm2','mean'),peri_density=('immune_density_peri_mm2','mean'),
        mean_boundary_coverage=('coliv_coverage_2um','mean'),mean_immune_localization=('immune_localization_fraction','mean')).reset_index();agg.to_csv(data/'all_h2_donor_core_peri_profiles.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(10.8,4.3),sharey=True)
    for ax,g in zip(axes,GROUPS):
        for r in agg[agg.donor_status.eq(g)].itertuples():ax.plot([0,1],[r.core_density,r.peri_density],'o-',color=COLORS[g],alpha=.7,lw=1,ms=4)
        ax.set(xticks=[0,1],xticklabels=['Core','20 µm peri'],title=g)
    axes[0].set_ylabel('Mean object immune density (cells/mm²)')
    save(fig,folder,'donor_core_peri_density', 'All2937 original eligible immune-exposed H2 objects are summarized within their20 donors (ND7, Aab+7, T1D6), including objects absent from the composition embedding. Each connected pair is the arithmetic mean of object-specific core/peri immune density for one donor. Counts, regional areas and object counts are retained in the source tables. The same cells may appear in multiple object peri associations; values are object means, not pooled unique tissue density. This is a donor profile with no paired test or additional hypothesis test.', ['all_h2_donor_core_peri_profiles.csv'])
    return dict(n_paired_objects=len(allq),n_embedded_complete_objects=len(q),n_donors=20,objects_by_disease=pop)


def main():
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    d,o=prepare();p9,support,r9=h9(d,o);p10,r10,lr=h10(d);clinical_plots(r10,lr);population=profiles(d,o)
    pd.DataFrame([p9,p10]).to_csv(OUT/'primary_comparisons.csv',index=False)
    original=(ROOT/'docs/T1D_research_hypotheses_20261022.md').read_text()
    conditions={}
    for h in ['H9','H10']:
        section=original.split('### '+h+' ',1)[1].split('\n### ',1)[0].split('\n## ',1)[0]
        conditions[h]={heading:section.split('**'+heading+'.** ',1)[1].split('\n\n',1)[0] for heading in ['If supported','If a meaningful difference is excluded'] if '**'+heading+'.** ' in section}
        if h=='H10':conditions[h]['If a meaningful improvement is excluded']=section.split('**If a meaningful improvement is excluded.** ',1)[1].split('\n\n',1)[0]
    (OUT/'original_if_then_conditions.json').write_text(json.dumps(conditions,indent=2)+'\n')
    audit=[dict(hypothesis='H9',item='Aab+ versus T1D H2-slope comparison adjusted for Beta quantity',execution=p9['status'],details=support,branch_decision='unassigned_user_judgment'),
        dict(hypothesis='H9',item='All13-donor conditional model and age/regional-immune sensitivities',execution='computed_supplementary_outside_shared_range',branch_decision='unassigned_user_judgment'),
        dict(hypothesis='H9',item='Single/multiple Aab, antibody specificity and continuous-age follow-ups',execution='age_adjusted_sensitivity_computed; antibody-specific comparisons not executed in this run',branch_decision='unassigned_user_judgment'),
        dict(hypothesis='H10',item='Original censored predictive log-score and conditional benchmark',execution='replaced_by_user_approved_nested_regression_and_held_out_RMSE',branch_decision='unassigned_user_judgment'),
        dict(hypothesis='H10',item='Nested regression, disease-conditioned sensitivity, proper censored likelihood, all six T1D',execution='computed',branch_decision='unassigned_user_judgment'),
        dict(hypothesis='H10',item='Original within-Beta INS tissue follow-up',execution='removed_by_user_instruction_no_quantitative_INS',branch_decision='unassigned_user_judgment')]
    (OUT/'if_then_execution_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (OUT/'figure_manifest.json').write_text(json.dumps(MANIFEST,indent=2)+'\n')
    pd.DataFrame(MANIFEST).to_csv(OUT/'figure_manifest.csv',index=False)
    primary=r10.iloc[0]
    text=['# Donor conditional and clinical measurements','',
        'Existing app clusters and pseudotime organize the donor profiles. They are not added to the fixed clinical model. All measurements, model terms and donor identities are exported; biological branch decisions remain unassigned.','',
        f'The coverage comparison uses the existing standardized H2 donor slopes. The common Beta-density range contains {p9["n_donors"]} donors and {p9["n_objects"]} eligible H2 objects; its adjusted design has {support["residual_df"]} residual degrees of freedom. Primary p value: unavailable. The separately reported all13-donor conditional coefficient is '+f'{r9.iloc[0].estimate:.6f} (95% CI {r9.iloc[0].lower_95:.6f} to {r9.iloc[0].upper_95:.6f}), p={r9.iloc[0].p_value:.6g}.','',
        f'The clinical primary has19 measured donors. Added-score F(1,{int(primary.df_denominator)})={primary.statistic:.6f}, p={primary.p_value:.6g}. Adjusted R²: baseline {primary.baseline_adjusted_r_squared:.6f}; augmented {primary.augmented_adjusted_r_squared:.6f}. Held-out RMSE: baseline {primary.baseline_held_out_rmse:.6f}; augmented {primary.augmented_held_out_rmse:.6f}, in natural-log C-peptide units. The F test and held-out error are distinct statistics.','',
        'The all20-donor sensitivity uses an actual censored-normal likelihood. The recorded <0.02 donor contributes a probability below the bound; it is never assigned a measured concentration or primary squared error. Disease-conditioned model results and all six T1D measurements/predictions are supplied.','',
        'Source-binding checks confirm that current clinical metadata, tissue-census cell/geometry inputs and H2 measurements are unchanged. H2 already rasterized the union of local capillaries before computing nonvascular coverage. The repaired geometric area sums and recovered structural-cell labels do not enter the fixed score.','',
        'The quantitative tissue-INS follow-up has been removed. Original if–then conditions are retained separately as historical text; no condition is declared activated.','',
        'Methods: approved nested ordinary regression/partial F test, donor-held-out RMSE, and properly censored-normal sensitivity; supplied-paper precedents include Drotar2025 PDFp7, Heide p23, Williams2024 p3, and supplied CELL-REPORTS draft p12. H9 reuses the original donor-slope conditional comparison, Barlow-style within-donor predictor scaling, and ordinary donor-level regression. No new hypothesis tests are applied to display-only cluster or pseudotime profiles.']
    (OUT/'README.md').write_text('\n'.join(text)+'\n')
    validation=dict(donors=20,groups=d.groupby('donor_status').size().to_dict(),primary_measured=19,censored_case=6550,all_six_T1D=True,
        no_quantitative_INS_endpoint=True,clinical_clusters_used_as_predictors=False,source_bindings_unchanged=True,profile_population=population,
        n_figures=len(MANIFEST),all_figures_have_three_formats=all((OUT/f).is_file() for m in MANIFEST for f in m['figure_files']))
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation),flush=True)


if __name__=='__main__':main()
