#!/usr/bin/env Rscript
# H4 predeclared model: raw Beta density, random donor intercept, continuous
# myeloid allocation, source signal and source/immune densities, size and class.
# Bodenmiller2025 PDF p28 supplies density LMMs; Drotar2025 p7 and Cellohood p23
# supply log/log1p transform precedents for the labelled transformation check.
suppressPackageStartupMessages(library(lmerTest))
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[[1]] else 'output/h1_h4_literature_methods_20261025/h4'
d <- read.csv(file.path(out,'object_source_features.csv'))
for(v in c('source_pool_nonzero','source_pool_complete','h1_primary'))
  d[[v]] <- tolower(as.character(d[[v]]))=='true'
d$case_id <- factor(d$case_id)
d$donor_status <- factor(d$donor_status,levels=c('ND','Aab+','T1D'))
d$object_class <- factor(d$object_class)
groups <- c('Pooled','ND','Aab+','T1D')
specs <- expand.grid(variant=c('continuous','lgals3_005','lgals3_010','lgals3_020',
                              'ck19_epithelial','all_exocrine'),cohort='all_complete',
                     outcome='raw_density',stringsAsFactors=FALSE)
specs <- rbind(specs,data.frame(variant='continuous',
     cohort=c('h1_primary','beta_positive','islets_only','complete_pool','all_complete'),
     outcome=c('raw_density','raw_density','raw_density','raw_density','log1p_density')))
write.csv(specs,file.path(out,'model_specifications_before_fit.csv'),row.names=FALSE)
results <- list(); details <- list(); allcoef <- list(); residuals <- list()
for(k in seq_len(nrow(specs))){
  spec<-specs[k,]; sub<-d[d$variant==spec$variant & d$source_pool_nonzero,]
  if(spec$cohort=='h1_primary')sub<-sub[sub$h1_primary,]
  if(spec$cohort=='beta_positive')sub<-sub[sub$n_beta_core>0,]
  if(spec$cohort=='islets_only')sub<-sub[sub$object_class=='Islet',]
  if(spec$cohort=='complete_pool')sub<-sub[sub$source_pool_complete,]
  for(g in groups){
    dd<-if(g=='Pooled')sub else sub[sub$donor_status==g,]
    dd<-droplevels(dd)
    dd$allocation_10pp<-dd$myeloid_fraction/.1
    dd$y<-if(spec$outcome=='raw_density')dd$beta_density_mm2 else log1p(dd$beta_density_mm2)
    covars<-c('signal_pool_density_mm2','myeloid_density_mm2','epithelial_density_mm2',
              'immune_density_mm2','core_area_um2')
    for(v in covars)dd[[paste0('z_',v)]]<-as.numeric(scale(dd[[v]]))
    terms<-c('allocation_10pp',paste0('z_',covars))
    if(nlevels(dd$object_class)>1)terms<-c(terms,'object_class')
    if(g=='Pooled')terms<-c(terms,'donor_status')
    form<-as.formula(paste('y ~',paste(terms,collapse=' + '),'+ (1|case_id)'))
    key<-paste(spec$variant,spec$cohort,spec$outcome,g,sep='__')
    warns<-character()
    fit<-tryCatch(withCallingHandlers(lmer(form,data=dd,REML=TRUE,
         control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),
         warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),
         error=function(e)e)
    if(inherits(fit,'error')){
      details[[key]]<-list(error=conditionMessage(fit));next
    }
    co<-coef(summary(fit));a<-co['allocation_10pp',]
    ci<-a['Estimate']+c(-1,1)*qt(.975,a['df'])*a['Std. Error']
    row<-data.frame(variant=spec$variant,cohort=spec$cohort,outcome=spec$outcome,
      group=g,objects=nrow(dd),donors=nlevels(dd$case_id),beta_positive=sum(dd$n_beta_core>0),
      beta_count=sum(dd$n_beta_core),estimate=a['Estimate'],se=a['Std. Error'],df=a['df'],
      ci_low=ci[1],ci_high=ci[2],p=a['Pr(>|t|)'],singular=isSingular(fit),
      converged=is.null(fit@optinfo$conv$lme4$messages),
      donor_sd=as.data.frame(VarCorr(fit))$sdcor[1],residual_sd=sigma(fit))
    results[[key]]<-row
    cdf<-as.data.frame(co);cdf$term<-rownames(co);cdf$model<-key;allcoef[[key]]<-cdf
    details[[key]]<-list(formula=deparse(form),warnings=warns,
         convergence=fit@optinfo$conv$lme4$messages,logLik=as.numeric(logLik(fit)))
    if(spec$variant=='continuous'&&spec$cohort=='all_complete'&&spec$outcome=='raw_density'){
      residuals[[key]]<-data.frame(islet_id=dd$islet_id,case_id=dd$case_id,group=g,
             observed=dd$y,fitted=fitted(fit),residual=residuals(fit),
             allocation_10pp=dd$allocation_10pp)
    }
    cat(key,' n=',nrow(dd),' estimate=',format(a['Estimate'],digits=5),
        ' p=',format(a['Pr(>|t|)'],digits=4),'\n',sep='')
  }
}
write.csv(do.call(rbind,results),file.path(out,'mixed_model_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,allcoef),file.path(out,'mixed_model_coefficients.csv'),row.names=FALSE)
write.csv(do.call(rbind,residuals),file.path(out,'primary_model_predictions.csv'),row.names=FALSE)
jsonlite::write_json(details,file.path(out,'mixed_model_details.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
