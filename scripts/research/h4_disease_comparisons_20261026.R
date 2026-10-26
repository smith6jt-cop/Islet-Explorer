#!/usr/bin/env Rscript
# Across-disease H4 interaction: Barlow PDF p22 random donor slopes/intercepts,
# within-donor standardized predictors and globally standardized outcome.
suppressPackageStartupMessages({library(lmerTest);library(emmeans);library(jsonlite)})
emm_options(lmerTest.limit=10000)
out <- 'output/h1_h4_disease_comparisons_20261026/h4'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source <- 'output/h1_h4_literature_methods_20261025/h4/object_source_features.csv'
raw <- read.csv(source,stringsAsFactors=FALSE)
for(v in c('source_pool_nonzero','source_pool_complete','h1_primary'))
  raw[[v]] <- tolower(as.character(raw[[v]]))=='true'
raw$case_id <- factor(raw$case_id)
raw$donor_status <- factor(raw$donor_status,levels=c('ND','Aab+','T1D'))
raw$object_class <- factor(raw$object_class)
specs<-data.frame(name=c('primary','lgals3_005','lgals3_010','lgals3_020',
                         'ck19_epithelial','all_exocrine','complete_pool','log1p_density'),
                  variant=c('continuous','lgals3_005','lgals3_010','lgals3_020',
                            'ck19_epithelial','all_exocrine','continuous','continuous'))
write.csv(specs,file.path(out,'specifications_before_fit.csv'),row.names=FALSE)
write_json(list(primary='Disease-by-source interaction and three BH-adjusted disease slope contrasts',
  predictor_standardization='Within donor z-scores for source fraction and each continuous adjustment',
  outcome_standardization='Beta density globally z-standardized across included objects',
  biological_replicate='Donor; random donor intercept and source slope',
  source=source,random_effects='(1 + source_z | case_id)',
  local_precedent='literature/Barlow2025.pdf PDF p22 and ~/Documents/papers/Bodenmiller2025.pdf PDF p28'),
  file.path(out,'analysis_specification.json'),pretty=TRUE,auto_unbox=TRUE)
covars<-c('signal_pool_density_mm2','myeloid_density_mm2','epithelial_density_mm2',
          'immune_density_mm2','core_area_um2')
all_slopes<-list();all_pairs<-list();all_omnibus<-list();diags<-list();scalers<-list()
for(k in seq_len(nrow(specs))){
  name<-specs$name[k];variant<-specs$variant[k]
  d<-droplevels(raw[raw$variant==variant & raw$source_pool_nonzero,])
  if(name=='complete_pool')d<-droplevels(d[d$source_pool_complete,])
  yy<-if(name=='log1p_density')log1p(d$beta_density_mm2) else d$beta_density_mm2
  d$beta_z<-as.numeric(scale(yy))
  for(v in c('myeloid_fraction',covars)){
    nm<-if(v=='myeloid_fraction')'source_z' else paste0(v,'_z')
    d[[nm]]<-ave(d[[v]],d$case_id,FUN=function(z){if(sd(z)==0)rep(0,length(z)) else as.numeric(scale(z))})
    scalers[[paste(name,v,sep='_')]]<-data.frame(model=name,variable=v,
       case_id=levels(d$case_id),mean=as.numeric(tapply(d[[v]],d$case_id,mean)),
       sd=as.numeric(tapply(d[[v]],d$case_id,sd)))
  }
  stopifnot(nlevels(d$case_id)==20,all(is.finite(d$beta_z)),
            all(tapply(d$source_z,d$case_id,sd)>.99))
  if(name=='primary'){
    stopifnot(nrow(d)==3681)
    write.csv(d,file.path(out,'primary_analysis_objects.csv'),row.names=FALSE)
  }
  form<-as.formula(paste('beta_z ~ donor_status * source_z +',
          paste(paste0(covars,'_z'),collapse=' + '),
          '+ object_class + (1 + source_z | case_id)'))
  warned<-character()
  fit<-withCallingHandlers(lmer(form,data=d,REML=TRUE,
          control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),
          warning=function(w){warned<<-c(warned,conditionMessage(w));invokeRestart('muffleWarning')})
  saveRDS(fit,file.path(out,paste0(name,'.rds')))
  trends<-emtrends(fit,~donor_status,var='source_z',lmer.df='satterthwaite')
  s<-as.data.frame(summary(trends,infer=c(TRUE,TRUE),adjust='none'))
  s$model<-name;s$donors<-as.numeric(table(d$donor_status)[s$donor_status])
  # objects and donors are reported separately; the latter are biological replicates.
  s$objects<-s$donors
  s$donors<-sapply(as.character(s$donor_status),function(g)nlevels(droplevels(d$case_id[d$donor_status==g])))
  all_slopes[[name]]<-s
  pair<-as.data.frame(summary(pairs(trends,adjust='none'),infer=c(TRUE,TRUE)))
  pair$p_bh<-p.adjust(pair$p.value,method='BH');pair$model<-name;all_pairs[[name]]<-pair
  aa<-as.data.frame(anova(fit,type=3,ddf='Satterthwaite'));aa$term<-rownames(aa);aa$model<-name
  all_omnibus[[name]]<-aa
  vc<-as.data.frame(VarCorr(fit));he<-eigen(fit@optinfo$derivs$Hessian,symmetric=TRUE,only.values=TRUE)$values
  diags[[name]]<-list(objects=nrow(d),donors=nlevels(d$case_id),formula=paste(deparse(form),collapse=' '),
     outcome_mean=mean(yy),outcome_sd=sd(yy),singular=isSingular(fit),
     optimizer_code=fit@optinfo$conv$opt,final_messages=fit@optinfo$conv$lme4$messages,
     warnings=warned,random_effects=vc,min_hessian_eigenvalue=min(he))
  cat('\n',name,' singular=',isSingular(fit),'\n',sep='');print(aa[grepl(':',aa$term),]);print(pair)
  if(name=='primary'){
     re<-coef(fit)$case_id
     coeff<-data.frame(case_id=rownames(re),donor_status=sapply(rownames(re),function(z)as.character(d$donor_status[d$case_id==z][1])),
                       slope=re$source_z)
     # coef() includes main source slope; add the disease-specific source interaction.
     fe<-fixef(fit)
     for(g in c('Aab+','T1D')){
       nm<-paste0('donor_status',g,':source_z')
       if(nm %in% names(fe)) coeff$slope[coeff$donor_status==g]<-coeff$slope[coeff$donor_status==g]+fe[nm]
     }
     write.csv(coeff,file.path(out,'donor_fitted_slopes.csv'),row.names=FALSE)
  }
}
write.csv(do.call(rbind,all_slopes),file.path(out,'disease_slopes.csv'),row.names=FALSE)
write.csv(do.call(rbind,all_pairs),file.path(out,'disease_slope_comparisons.csv'),row.names=FALSE)
write.csv(do.call(rbind,all_omnibus),file.path(out,'model_anova.csv'),row.names=FALSE)
write.csv(do.call(rbind,scalers),file.path(out,'within_donor_scaling.csv'),row.names=FALSE)
write_json(diags,file.path(out,'model_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
