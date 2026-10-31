#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(lmerTest);library(emmeans);library(jsonlite)})
emm_options(lmerTest.limit=10000)
out<-'output/h5_h10_original_plan_20261026/h5'
raw<-read.csv(file.path(out,'object_ins_features.csv'),stringsAsFactors=FALSE)
raw$case_id<-factor(raw$case_id);raw$donor_status<-factor(raw$donor_status,levels=c('ND','Aab+','T1D'))
raw$object_class<-factor(raw$object_class)
raw$h5_primary_eligible<-tolower(as.character(raw$h5_primary_eligible))=='true'
raw$log_beta_density<-log1p(raw$beta_density_mm2)
raw$log_core_area<-log(raw$core_area_um2)
raw$log_nearby<-log1p(raw$n_immune_nearby)
spec<-data.frame(model=c('primary','median_INS','raw_INS','corrected_INS','strong_hormone_margin'),
  outcome=c('ins_mean_beta_core','ins_median_beta_core','ins_mean_beta_raw',
            'ins_mean_beta_corrected','ins_mean_beta_margin4'))
write.csv(spec,file.path(out,'model_specifications_before_fit.csv'),row.names=FALSE)
write_json(list(primary='Disease comparison of donor-specific immune-penetration slopes for within-Beta mean INS; donor regression used after singular random-slope fits',
  primary_adjustment=c('Beta density','core area','object class','donor/image intercept'),
  primary_formula='within-Beta mean INS z ~ immune-penetration z + log1p(Beta density) z + log(core area) z + object class; fitted separately for each donor',
  donor_regression_sensitivity='Adds log1p(nearby immune burden) z, retaining all nonconstant terms; full model rank and positive residual degrees of freedom required',
  mixed_model_role='Separately adjusted supporting analyses, including nearby immune burden in every fit; not identical to primary donor-regression covariates',
  mixed_model_adjustment=c('Beta density','core area','nearby immune burden','object class','donor/image intercept'),
  background='One source image per donor; donor-specific intercept accounts for image-constant background. Measured epithelial-negative background values retained in donor_background.csv.',
  standardization='Predictors within donor, outcome across included objects',
  no_beta_handling='Within-Beta mean is unavailable, never zero',
  correlated_random_slopes='(1+exposure_z|case_id)',
  singular_sensitivity='If correlated model is singular, independent intercept/slope sensitivity is retained separately'),
  file.path(out,'analysis_specification.json'),pretty=TRUE,auto_unbox=TRUE)
slopes<-list();pairs<-list();tests<-list();diags<-list()
for(k in seq_len(nrow(spec))){
  name<-spec$model[k];outcome<-spec$outcome[k]
  d<-droplevels(raw[raw$h5_primary_eligible&is.finite(raw[[outcome]]),])
  d$y_z<-as.numeric(scale(d[[outcome]]))
  for(v in c('immune_localization_fraction','log_beta_density','log_core_area','log_nearby')){
    nm<-if(v=='immune_localization_fraction')'exposure_z' else paste0(v,'_z')
    d[[nm]]<-ave(d[[v]],d$case_id,FUN=function(z)if(sd(z)==0)rep(0,length(z)) else as.numeric(scale(z)))
  }
  if(name=='primary')write.csv(d,file.path(out,'primary_analysis_objects.csv'),row.names=FALSE)
  fixed<-'y_z ~ exposure_z*donor_status + log_beta_density_z + log_core_area_z + log_nearby_z + object_class'
  randoms<-c(correlated='(1+exposure_z|case_id)')
  for(random_name in c('correlated','independent')){
    if(random_name=='independent'&&!singular)next
    rand<-if(random_name=='correlated')'(1+exposure_z|case_id)' else '(1|case_id)+(0+exposure_z|case_id)'
    form<-as.formula(paste(fixed,'+',rand));key<-paste(name,random_name,sep='_')
    fit<-lmer(form,data=d,REML=TRUE,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000)))
    singular<-isSingular(fit)
    saveRDS(fit,file.path(out,paste0(key,'.rds')))
    aa<-as.data.frame(anova(fit,type=3,ddf='Satterthwaite'));aa$term<-rownames(aa);aa$model<-key;tests[[key]]<-aa
    tr<-emtrends(fit,~donor_status,var='exposure_z',lmer.df='satterthwaite')
    ss<-as.data.frame(summary(tr,infer=c(TRUE,TRUE),adjust='none'));ss$model<-key
    ss$objects<-as.numeric(table(d$donor_status)[ss$donor_status]);ss$donors<-sapply(as.character(ss$donor_status),function(g)nlevels(droplevels(d$case_id[d$donor_status==g])))
    slopes[[key]]<-ss
    pp<-as.data.frame(summary(pairs(tr,adjust='none'),infer=c(TRUE,TRUE)));pp$p_bh<-p.adjust(pp$p.value,'BH');pp$model<-key;pairs[[key]]<-pp
    vc<-as.data.frame(VarCorr(fit));he<-eigen(fit@optinfo$derivs$Hessian,symmetric=TRUE,only.values=TRUE)$values
    diags[[key]]<-list(objects=nrow(d),donors=nlevels(d$case_id),outcome=outcome,
       singular=singular,optimizer_code=fit@optinfo$conv$opt,final_messages=fit@optinfo$conv$lme4$messages,
       random_covariance=vc,min_hessian_eigenvalue=min(he),formula=paste(deparse(form),collapse=' '))
    cat('\n',key,' singular=',singular,'\n',sep='');print(aa[grepl(':',aa$term),]);print(pp)
  }
}
write.csv(do.call(rbind,slopes),file.path(out,'disease_slopes.csv'),row.names=FALSE)
write.csv(do.call(rbind,pairs),file.path(out,'disease_slope_comparisons.csv'),row.names=FALSE)
write.csv(do.call(rbind,tests),file.path(out,'model_anova.csv'),row.names=FALSE)
write_json(diags,file.path(out,'model_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
