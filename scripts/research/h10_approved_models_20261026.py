#!/usr/bin/env python3
"""Approved H10: nested donor regressions, held-out error, censored sensitivity.

The fixed score and predictors are specified before fitting. Clinical censoring
is retained; no invented concentration replaces the below-limit observation.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import optimize, stats
import statsmodels.api as sm
from statsmodels.tools.numdiff import approx_hess

REPO=Path(__file__).resolve().parents[2]
OUT=REPO/'output/h5_h10_original_plan_20261026/h10'
SOURCE=OUT/'donor_inputs.csv'
H5=REPO/'output/h5_h10_original_plan_20261026/h5/plot_data_donor_ins_measurements.csv'
BASE=['tissue_beta_density_mm2','age_years_npod','admission_glucose_mg_dL']
COMPONENTS=['mean_boundary_coverage','mean_immune_localization']
GROUPS=['ND','Aab+','T1D']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def design(training,target,augmented,disease):
    """Every learned scale is calculated on training donors only."""
    columns=BASE+COMPONENTS
    mu=training[columns].mean(); sd=training[columns].std(ddof=1)
    assert np.isfinite(mu).all() and (sd>0).all()
    scaled=(target[columns]-mu)/sd
    x=pd.DataFrame({'Intercept':np.ones(len(target))},index=target.index)
    for col in BASE:x[col+'_z']=scaled[col]
    if disease:
        for group in GROUPS[1:]:x['disease_'+group]=target.donor_status.eq(group).astype(float)
    if augmented:x['spatial_score']=.5*(scaled.mean_boundary_coverage-scaled.mean_immune_localization)
    scaling=pd.DataFrame({'variable':columns,'training_mean':mu.values,'training_sd':sd.values})
    return x,scaling


def ordinary_analysis(data,outcome,label,disease=False):
    prefix=label+('_disease_conditioned' if disease else '')
    fits={};metrics=[];coefs=[];full=[]
    for augmented in [False,True]:
        model='augmented' if augmented else 'baseline'
        x,scaling=design(data,data,augmented,disease)
        assert np.linalg.matrix_rank(x)==x.shape[1]
        m=sm.OLS(data[outcome],x).fit()
        fits[model]=m;m.save(OUT/f'{prefix}_{model}.pickle')
        ci=m.conf_int()
        for term in x.columns:
            coefs.append(dict(analysis=prefix,model=model,term=term,estimate=m.params[term],SE=m.bse[term],
                lower_95=ci.loc[term,0],upper_95=ci.loc[term,1],p_value=m.pvalues[term],residual_df=m.df_resid))
        metrics.append(dict(analysis=prefix,model=model,donors=len(data),parameters=x.shape[1],
            residual_df=m.df_resid,r_squared=m.rsquared,adjusted_r_squared=m.rsquared_adj,
            in_sample_rmse=float(np.sqrt(np.mean(m.resid**2))),residual_sigma=float(np.sqrt(m.mse_resid))))
        q=data[['case_id','donor_status',outcome]].copy();q['model']=model;q['fitted']=m.fittedvalues;q['residual']=m.resid
        q['leverage']=m.get_influence().hat_matrix_diag
        if augmented:q['spatial_score']=x.spatial_score
        full.append(q)
    pd.DataFrame(coefs).to_csv(OUT/f'{prefix}_coefficients.csv',index=False)
    pd.concat(full).to_csv(OUT/f'{prefix}_fitted_values.csv',index=False)
    scaling.to_csv(OUT/f'{prefix}_full_sample_scaling.csv',index=False)
    f,p,df_num=fits['augmented'].compare_f_test(fits['baseline'])

    predictions=[];foldscales=[];foldcoefs=[];foldchecks=[]
    for case in data.case_id:
        training=data.loc[data.case_id.ne(case)];test=data.loc[data.case_id.eq(case)]
        row=dict(case_id=int(case),donor_status=test.donor_status.iloc[0],observed=float(test[outcome].iloc[0]),
                 analysis=prefix,training_donors=len(training))
        for augmented in [False,True]:
            model='augmented' if augmented else 'baseline'
            x,sc=design(training,training,augmented,disease)
            xt,_=design(training,test,augmented,disease)
            assert np.linalg.matrix_rank(x)==x.shape[1]
            fit=sm.OLS(training[outcome],x).fit()
            pred=float(fit.predict(xt).iloc[0]);row[model+'_prediction']=pred
            row[model+'_error']=pred-row['observed'];row[model+'_squared_error']=row[model+'_error']**2
            if augmented:
                row['spatial_score_training_scaled']=float(xt.spatial_score.iloc[0])
                sc['held_out_case_id']=int(case);sc['training_donors']=len(training);foldscales.append(sc)
            for term,value in fit.params.items():foldcoefs.append(dict(held_out_case_id=case,model=model,term=term,coefficient=value))
            foldchecks.append(dict(held_out_case_id=case,model=model,training_donors=len(training),
                held_out_absent_from_training=case not in training.case_id.values,full_rank=True,residual_df=fit.df_resid))
        row['augmented_minus_baseline_squared_error']=row['augmented_squared_error']-row['baseline_squared_error']
        predictions.append(row)
    pred=pd.DataFrame(predictions)
    pred.to_csv(OUT/f'{prefix}_held_out_predictions.csv',index=False)
    pd.concat(foldscales).to_csv(OUT/f'{prefix}_held_out_scaling.csv',index=False)
    pd.DataFrame(foldcoefs).to_csv(OUT/f'{prefix}_held_out_coefficients.csv',index=False)
    pd.DataFrame(foldchecks).to_csv(OUT/f'{prefix}_held_out_validation.csv',index=False)
    metrics=pd.DataFrame(metrics)
    for model in ['baseline','augmented']:
        metrics.loc[metrics.model.eq(model),'held_out_rmse']=np.sqrt(pred[model+'_squared_error'].mean())
        metrics.loc[metrics.model.eq(model),'held_out_mae']=pred[model+'_error'].abs().mean()
    metrics.to_csv(OUT/f'{prefix}_model_metrics.csv',index=False)
    groups=[]
    for group,z in [('All',pred)]+list(pred.groupby('donor_status')):
        for model in ['baseline','augmented']:
            groups.append(dict(analysis=prefix,group=group,model=model,donors=len(z),rmse=np.sqrt(z[model+'_squared_error'].mean()),
                               mae=z[model+'_error'].abs().mean()))
    pd.DataFrame(groups).to_csv(OUT/f'{prefix}_held_out_error_by_group.csv',index=False)
    bm=metrics.set_index('model').loc['baseline'];am=metrics.set_index('model').loc['augmented']
    comparison=dict(analysis=prefix,donors=len(data),test='Partial F test for added fixed spatial score',
        statistic=float(f),df_numerator=float(df_num),df_denominator=float(fits['augmented'].df_resid),p_value=float(p),
        score_coefficient=float(fits['augmented'].params.spatial_score),
        score_lower_95=float(fits['augmented'].conf_int().loc['spatial_score',0]),
        score_upper_95=float(fits['augmented'].conf_int().loc['spatial_score',1]),
        baseline_adjusted_r_squared=float(bm.adjusted_r_squared),augmented_adjusted_r_squared=float(am.adjusted_r_squared),
        baseline_held_out_rmse=float(bm.held_out_rmse),augmented_held_out_rmse=float(am.held_out_rmse),
        augmented_minus_baseline_rmse=float(am.held_out_rmse-bm.held_out_rmse))
    print(json.dumps(comparison),flush=True)
    return comparison,pred,fits


def censored_fit(data,augmented,disease):
    x,scaling=design(data,data,augmented,disease);names=list(x);X=x.to_numpy()
    y=np.log(data.C_peptide_ng_mL.to_numpy());c=data.C_peptide_censored.to_numpy(bool)
    assert np.linalg.matrix_rank(X)==X.shape[1] and c.sum()==1
    # The censored row contributes P(log C-peptide < log 0.02), not a density at the bound.
    def value_gradient(theta):
        beta=theta[:-1];sigma=np.exp(theta[-1]);r=(y-X@beta)/sigma
        ll=np.where(c,stats.norm.logcdf(r),stats.norm.logpdf(r)-np.log(sigma))
        mills=np.exp(stats.norm.logpdf(r[c])-stats.norm.logcdf(r[c]))
        weight=-r/sigma;weight[c]=mills/sigma
        glog=1-r*r;glog[c]=mills*r[c]
        return -float(ll.sum()),np.r_[X.T@weight,glog.sum()]
    initial=sm.OLS(y[~c],X[~c]).fit()
    theta0=np.r_[initial.params,np.log(np.sqrt(np.mean(initial.resid**2)))]
    fit=optimize.minimize(value_gradient,theta0,jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':10000})
    # Retry numerical convergence only; both candidates optimize the same fixed likelihood.
    if np.max(np.abs(value_gradient(fit.x)[1]))>1e-6:
        fit=optimize.minimize(value_gradient,fit.x,jac=True,method='L-BFGS-B',
            bounds=[(None,None)]*(len(theta0)-1)+[(-10,10)],options={'ftol':1e-14,'gtol':1e-9,'maxiter':10000})
    nll,gradient=value_gradient(fit.x)
    hessian=approx_hess(fit.x,lambda a:value_gradient(a)[0]);eig=np.linalg.eigvalsh(hessian)
    assert np.max(np.abs(gradient))<1e-5 and eig.min()>0
    covariance=np.linalg.inv(hessian);se=np.sqrt(np.diag(covariance));sigma=np.exp(fit.x[-1])
    model='augmented' if augmented else 'baseline';label='censored'+('_disease_conditioned' if disease else '')
    rows=[]
    for j,term in enumerate(names+['log_residual_sigma']):
        rows.append(dict(analysis=label,model=model,term=term,estimate=fit.x[j],SE=se[j],
            lower_95=fit.x[j]-stats.norm.ppf(.975)*se[j],upper_95=fit.x[j]+stats.norm.ppf(.975)*se[j],
            p_value=2*stats.norm.sf(abs(fit.x[j]/se[j]))))
    predictions=data[['case_id','donor_status','C_peptide_ng_mL','C_peptide_censored']].copy()
    predictions['analysis']=label;predictions['model']=model;predictions['predicted_log_c_peptide']=X@fit.x[:-1]
    predictions['predicted_median_c_peptide']=np.exp(predictions.predicted_log_c_peptide)
    predictions['probability_below_0_02']=stats.norm.cdf((np.log(.02)-predictions.predicted_log_c_peptide)/sigma)
    if augmented:predictions['spatial_score']=x.spatial_score.values
    diagnostic=dict(analysis=label,model=model,donors=len(data),uncensored_donors=int((~c).sum()),censored_donors=int(c.sum()),
        log_likelihood=-nll,parameters=len(fit.x),aic=2*len(fit.x)+2*nll,residual_sigma=sigma,
        optimizer_success=bool(fit.success),optimizer_message=str(fit.message),
        maximum_absolute_gradient=float(np.max(np.abs(gradient))),minimum_hessian_eigenvalue=float(eig.min()),full_rank=True)
    (OUT/f'{label}_{model}_fit.json').write_text(json.dumps(dict(diagnostics=diagnostic,terms=names+['log_residual_sigma'],
        parameters=fit.x.tolist(),covariance=covariance.tolist(),likelihood='normal log density for measured values; normal log CDF at recorded upper bound for censored values'),indent=2)+'\n')
    scaling.to_csv(OUT/f'{label}_{model}_scaling.csv',index=False)
    return rows,predictions,diagnostic


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source_hash=sha(SOURCE);h5_hash=sha(H5)
    d=pd.read_csv(SOURCE).sort_values('case_id').reset_index(drop=True)
    assert len(d)==20 and d.case_id.is_unique and 78 not in d.case_id.values
    assert d[BASE+COMPONENTS+['C_peptide_ng_mL']].notna().all().all()
    assert d.C_peptide_censored.sum()==1 and d.loc[d.C_peptide_censored,'case_id'].tolist()==[6550]
    d['log_c_peptide']=np.log(d.C_peptide_ng_mL)
    measured=d[~d.C_peptide_censored].copy()
    specification=dict(status='User-approved conventional analysis with proper censored sensitivity',
        input_files={str(SOURCE):source_hash,str(H5):h5_hash},
        primary_outcome='Natural log measured C-peptide (19 donors); no pseudocount',baseline=BASE,
        score='0.5*(z mean_boundary_coverage - z mean_immune_localization); equal weights fixed before clinical fitting',
        primary_test='Nested ordinary regression partial F test for added score; adjusted R-squared',
        prediction='Leave one donor out; training-only scales; RMSE on natural-log C-peptide; no imputation or tuning',
        disease_sensitivity='Add Aab+ and T1D indicators to baseline and augmented models',
        censoring_sensitivity='All20 donors: normal likelihood for log C-peptide with log CDF for <0.02; baseline-versus-augmented likelihood ratio',
        tissue_followup='H5 existing donor cell-weighted mean INS within typed Beta cells; same fixed baseline and score; raw scaled INS outcome with held-out RMSE',
        T1D_display='All six donors: five measured outcomes with primary LOO predictions; censored6550 predicted from all19 measured donors and never assigned a squared error',
        inferential_scope='Donors are replicates. The primary association test and held-out error are separate statistics. Censored fits and tissue follow-up are secondary. Branch decision: unassigned (user judgment).')
    (OUT/'approved_analysis_specification.json').write_text(json.dumps(specification,indent=2)+'\n')
    comparisons=[]
    primary,pred,fits=ordinary_analysis(measured,'log_c_peptide','clinical',False);comparisons.append(primary)
    disease,_,_=ordinary_analysis(measured,'log_c_peptide','clinical',True);comparisons.append(disease)

    # One genuinely withheld censored donor: predictions use all19 measured donors, without imputing its outcome.
    censored=d[d.C_peptide_censored].copy();censor_prediction=censored[['case_id','donor_status','C_peptide_ng_mL','C_peptide_censored']].copy()
    for aug in [False,True]:
        model='augmented' if aug else 'baseline';xt,_=design(measured,censored,aug,False)
        censor_prediction[model+'_prediction']=fits[model].predict(xt).to_numpy()
    censor_prediction['prediction_basis']='Model trained on all19 measured donors; censored donor was never in training'
    censor_prediction.to_csv(OUT/'censored_donor_unscored_primary_predictions.csv',index=False)
    t1d=d[d.donor_status.eq('T1D')][['case_id','donor_status','C_peptide_ng_mL','C_peptide_censored']].merge(
        pd.concat([pred[['case_id','baseline_prediction','augmented_prediction']],censor_prediction[['case_id','baseline_prediction','augmented_prediction']]]),
        on='case_id',validate='one_to_one')
    for model in ['baseline','augmented']:
        t1d[model+'_predicted_median_ng_mL']=np.exp(t1d[model+'_prediction'])
        t1d[model+'_log_error']=np.where(t1d.C_peptide_censored,np.nan,t1d[model+'_prediction']-np.log(t1d.C_peptide_ng_mL))
    t1d['prediction_basis']=np.where(t1d.C_peptide_censored,'Trained on19 measured donors; outcome unscored','Leave-one-donor-out; trained on18 measured donors')
    t1d.to_csv(OUT/'all_six_T1D_observations_and_predictions.csv',index=False)

    cens_rows=[];cens_predictions=[];diagnostics=[];likelihood_tests=[]
    for conditioned in [False,True]:
        info=[]
        for augmented in [False,True]:
            rows,values,diagnostic=censored_fit(d,augmented,conditioned)
            cens_rows+=rows;cens_predictions.append(values);diagnostics.append(diagnostic);info.append(diagnostic)
        lr=2*(info[1]['log_likelihood']-info[0]['log_likelihood'])
        assert lr>=-1e-8
        likelihood_tests.append(dict(analysis=info[0]['analysis'],donors=20,censored_donors=1,
            test='Likelihood-ratio test for added score in censored-normal log-C-peptide regression',
            statistic=max(0.,lr),df=1,p_value=stats.chi2.sf(max(0.,lr),1),
            baseline_log_likelihood=info[0]['log_likelihood'],augmented_log_likelihood=info[1]['log_likelihood']))
    pd.DataFrame(cens_rows).to_csv(OUT/'censored_model_coefficients.csv',index=False)
    pd.concat(cens_predictions).to_csv(OUT/'censored_model_fitted_values.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(OUT/'censored_model_diagnostics.csv',index=False)
    pd.DataFrame(likelihood_tests).to_csv(OUT/'censored_likelihood_comparisons.csv',index=False)

    h5=pd.read_csv(H5)
    tissue=d.merge(h5[['case_id','mean_beta_INS','n_beta_ins_measured','beta_positive_objects']],on='case_id',validate='one_to_one')
    assert tissue.mean_beta_INS.notna().all() and tissue.n_beta_ins_measured.gt(0).all()
    tissue.to_csv(OUT/'tissue_followup_inputs.csv',index=False)
    follow,_,_=ordinary_analysis(tissue,'mean_beta_INS','tissue_INS',False);comparisons.append(follow)
    follow_dis,_,_=ordinary_analysis(tissue,'mean_beta_INS','tissue_INS',True);comparisons.append(follow_dis)
    pd.DataFrame(comparisons).to_csv(OUT/'ordinary_model_comparisons.csv',index=False)
    pd.DataFrame([primary]).to_csv(OUT/'primary_test.csv',index=False)
    validation=dict(donors=20,primary_measured_donors=19,primary_T1D_donors=5,censored_case_id=6550,
        held_out_preprocessing_training_only=True,imputed_values=0,censor_bound_substitutions=0,
        censored_model_max_gradient=max(x['maximum_absolute_gradient'] for x in diagnostics),
        censored_model_min_hessian_eigenvalue=min(x['minimum_hessian_eigenvalue'] for x in diagnostics),
        source_unchanged=sha(SOURCE)==source_hash,H5_source_unchanged=sha(H5)==h5_hash,
        script_sha256=sha(__file__),source_sha256=source_hash,H5_source_sha256=h5_hash)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(pd.DataFrame(likelihood_tests).to_string(index=False),flush=True)
    print(json.dumps(validation),flush=True)


if __name__=='__main__':main()
