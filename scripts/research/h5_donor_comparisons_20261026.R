#!/usr/bin/env Rscript
# Donor-unit conditional comparison after the primary random slope covariance
# reached its boundary. Preserve the H5 density/size/class adjustments.
out<-'output/h5_h10_original_plan_20261026/h5'
d<-read.csv(file.path(out,'primary_analysis_objects.csv'))
d$case_id<-factor(d$case_id);d$donor_status<-factor(d$donor_status,levels=c('ND','Aab+','T1D'))
d$object_class<-factor(d$object_class)
for(variant in c('original_covariates','additional_nearby_burden')){
prefix<-if(variant=='original_covariates')'' else 'nearby_burden_sensitivity_'
rows<-list();fits<-list()
for(case in levels(d$case_id)){
  z<-droplevels(d[d$case_id==case,])
  controls<-c('log_beta_density_z','log_core_area_z')
  if(variant=='additional_nearby_burden')controls<-c(controls,'log_nearby_z')
  controls<-controls[vapply(controls,function(v)sd(z[[v]])>0,logical(1))]
  if(nlevels(z$object_class)>1)controls<-c(controls,'object_class')
  form<-as.formula(paste('y_z ~ exposure_z +',paste(controls,collapse=' + ')))
  fit<-lm(form,data=z)
  estimable<-fit$rank==ncol(model.matrix(fit))&&df.residual(fit)>0&&is.finite(coef(fit)['exposure_z'])
  rr<-data.frame(case_id=case,donor_status=as.character(z$donor_status[1]),objects=nrow(z),
     model_parameters=ncol(model.matrix(fit)),model_rank=fit$rank,residual_df=df.residual(fit),
     estimable=estimable,standardized_beta=NA_real_,SE=NA_real_,lower=NA_real_,upper=NA_real_,
     formula=paste(deparse(form),collapse=' '))
  if(estimable){
    co<-coef(summary(fit))['exposure_z',];ci<-confint(fit,'exposure_z')
    rr$standardized_beta<-co['Estimate'];rr$SE<-co['Std. Error'];rr$lower<-ci[1];rr$upper<-ci[2]
  }
  fits[[case]]<-fit;rows[[case]]<-rr
}
x<-do.call(rbind,rows)
write.csv(x,file.path(out,paste0(prefix,'donor_regression_eligibility.csv')),row.names=FALSE)
saveRDS(fits,file.path(out,paste0(prefix,'donor_linear_models.rds')))
x<-x[x$estimable,];write.csv(x,file.path(out,paste0(prefix,'plot_data_donor_coefficients.csv')),row.names=FALSE)
kw<-kruskal.test(standardized_beta~donor_status,data=x)
write.csv(data.frame(test='Kruskal-Wallis',statistic=unname(kw$statistic),df=unname(kw$parameter),
  p_value=kw$p.value,donors=nrow(x)),file.path(out,paste0(prefix,'donor_disease_global_test.csv')),row.names=FALSE)
tests<-list()
for(p in combn(c('ND','Aab+','T1D'),2,simplify=FALSE)){
  a<-x$standardized_beta[x$donor_status==p[1]];b<-x$standardized_beta[x$donor_status==p[2]]
  tied<-anyDuplicated(c(a,b))>0
  w<-wilcox.test(a,b,exact=!tied,correct=TRUE)
  tests[[paste(p,collapse='_')]]<-data.frame(group_1=p[1],group_2=p[2],donors_1=length(a),donors_2=length(b),
      p_value=w$p.value,method=if(tied)'Tie-corrected Wilcoxon rank-sum' else 'Exact Wilcoxon rank-sum')
}
t<-do.call(rbind,tests);t$p_bh<-p.adjust(t$p_value,'BH')
write.csv(t,file.path(out,paste0(prefix,'donor_disease_pairwise_tests.csv')),row.names=FALSE)
cat('\nVariant:',variant,'\n')
print(kw);print(t);print(do.call(rbind,rows)[!do.call(rbind,rows)$estimable,])
}
