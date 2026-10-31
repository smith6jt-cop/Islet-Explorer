#!/usr/bin/env Rscript
# New non-INS H7 outcome; protocol_before_fit.json precedes all fits.
options(warn=1)
suppressPackageStartupMessages({library(lmerTest);library(emmeans);library(jsonlite)})
emm_options(lmerTest.limit=10000)
out<-'output/app_hypotheses_20261026/h6_h8';models<-file.path(out,'models')
stopifnot(file.exists(file.path(out,'protocol_before_fit.json')))
d<-read.csv(file.path(out,'data/h7_objects.csv'),stringsAsFactors=FALSE)
d$case_id<-factor(d$case_id);d$donor_status<-factor(d$donor_status,levels=c('ND','Aab+','T1D'))
truth<-function(x)tolower(as.character(x))=='true'
controls<-c('nerve_area_fraction','cd31_area_fraction','log_core_area','log_perimeter','log_nearby_immune_density')
spec<-data.frame(model=c('primary_20um','distance_10um','distance_5um','overlap','GAP43_20um'),
 exposure=c('nerve_near_cd31_20um','nerve_near_cd31_10um','nerve_near_cd31_5um','nerve_vessel_overlap','gap43_nerve_near_cd31_20um'),
 outcome='log1p_beta_density',role=c('primary',rep('sensitivity',4)))
write.csv(spec,file.path(models,'h7_model_specifications.csv'),row.names=FALSE)
within_z<-function(x,g)ave(x,g,FUN=function(v)if(length(v)>1&&sd(v)>1e-12)(v-mean(v))/sd(v)else rep(0,length(v)))
bind<-function(x)if(length(x))do.call(rbind,x)else data.frame()
slopes<-list();pairs_out<-list();global<-list();donor_rows<-list();diag<-list();eligibility<-list();fixed<-list();scaling<-list()
for(i in seq_len(nrow(spec))){
 name<-spec$model[i];exposure<-spec$exposure[i]
 use<-truth(d$core_complete)&truth(d$neighborhood_complete_20)&is.finite(d[[exposure]])&is.finite(d$log1p_beta_density)
 for(v in controls)use<-use&is.finite(d[[v]])
 z<-droplevels(d[use,]);z$exposure_z<-within_z(z[[exposure]],z$case_id)
 z$response_z<-as.numeric(scale(z$log1p_beta_density))
 for(v in controls)z[[paste0(v,'_z')]]<-within_z(z[[v]],z$case_id)
 control_names<-paste0(controls,'_z')
 write.csv(z,file.path(models,paste0('h7_',name,'_objects.csv')),row.names=FALSE)
 scaling[[name]]<-data.frame(model=name,outcome_mean=mean(z$log1p_beta_density),outcome_sd=sd(z$log1p_beta_density),n_objects=nrow(z),n_donors=nlevels(z$case_id),zero_Beta_objects=sum(z$n_beta_core==0))
 eligibility[[name]]<-bind(lapply(levels(d$case_id),function(ca){
  all<-d[d$case_id==ca,];zz<-z[z$case_id==ca,]
  data.frame(model=name,case_id=ca,donor_status=as.character(all$donor_status[1]),all_objects=nrow(all),model_objects=nrow(zz),zero_Beta_objects=sum(zz$n_beta_core==0),exposure_sd=if(nrow(zz)>1)sd(zz[[exposure]])else NA_real_)}))
 current_donors<-list();donor_models<-list()
 for(ca in levels(d$case_id)){
  zz<-z[z$case_id==ca,];status<-as.character(d$donor_status[d$case_id==ca][1])
  form<-reformulate(c('exposure_z',control_names),response='response_z')
  fm<-if(nrow(zz))lm(form,zz)else NULL
  varying<-nrow(zz)>1&&sd(zz$exposure_z)>1e-12
  valid<-!is.null(fm)&&varying&&fm$rank==ncol(model.matrix(fm))&&df.residual(fm)>0&&is.finite(coef(fm)['exposure_z'])
  constant_response<-nrow(zz)>1&&sd(zz$response_z)<1e-12
  cf<-if(valid)coef(summary(fm))['exposure_z',]else rep(NA_real_,4)
  ci<-if(valid)confint(fm,'exposure_z')else matrix(c(NA_real_,NA_real_),nrow=1)
  if(valid&&constant_response){cf[1:2]<-0;ci[,]<-0}
  current_donors[[ca]]<-data.frame(model=name,case_id=ca,donor_status=status,objects=nrow(zz),estimable=valid,
   exposure_varies=varying,constant_response=constant_response,design_rank=if(is.null(fm))0 else fm$rank,
   design_columns=if(is.null(fm))NA else ncol(model.matrix(fm)),residual_df=if(is.null(fm))NA else df.residual(fm),
   estimate=cf[1],SE=cf[2],lower=ci[1],upper=ci[2],
   reason=if(!nrow(zz))'no complete measurements'else if(!varying)'constant exposure'else if(fm$rank<ncol(model.matrix(fm)))'full covariate design not estimable'else if(df.residual(fm)<=0)'no residual degrees of freedom'else if(constant_response)'constant outcome; exact zero slope'else'estimable')
  donor_models[[ca]]<-fm
 }
 dc<-bind(current_donors);donor_rows[[name]]<-dc;saveRDS(donor_models,file.path(models,paste0('h7_',name,'_donor_models.rds')))
 ff<-paste('response_z ~ exposure_z * donor_status +',paste(control_names,collapse=' + '))
 chosen<-NULL;attempts<-list()
 for(covariance in c('correlated','independent')){
  random<-if(covariance=='correlated')'(1+exposure_z|case_id)'else'(1+exposure_z||case_id)'
  warning_text<-character();form<-as.formula(paste(ff,'+',random))
  fit<-tryCatch(withCallingHandlers(lmer(form,z,REML=TRUE,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),warning=function(w){warning_text<<-c(warning_text,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error')){attempts[[covariance]]<-list(error=conditionMessage(fit));next}
  saveRDS(fit,file.path(models,paste0('h7_',name,'_',covariance,'.rds')))
  singular<-isSingular(fit);conv<-fit@optinfo$conv$lme4$messages
  attempts[[covariance]]<-list(singular=singular,optimizer_code=fit@optinfo$conv$opt,convergence_messages=conv,warnings=warning_text,variance=as.data.frame(VarCorr(fit)),formula=paste(deparse(form),collapse=' '))
  if(!singular&&length(conv)==0&&fit@optinfo$conv$opt==0){chosen<-fit;chosen_covariance<-covariance;break}
 }
 diag[[name]]<-attempts
 extract_fit<-function(fit,label,role){
  tr<-emtrends(fit,~donor_status,var='exposure_z',lmer.df='satterthwaite')
  ss<-as.data.frame(summary(tr,infer=c(TRUE,TRUE),adjust='none'))
  sr<-data.frame(model=name,source=label,role=role,donor_status=as.character(ss$donor_status),estimate=ss$exposure_z.trend,SE=ss$SE,df=ss$df,lower=ss$lower.CL,upper=ss$upper.CL,p_value=ss$p.value,p_bh=p.adjust(ss$p.value,'BH',n=3),n_donors=nlevels(z$case_id),n_objects=nrow(z))
  pp<-as.data.frame(summary(pairs(tr),infer=c(TRUE,TRUE),adjust='none'))
  pr<-data.frame(model=name,source=label,role=role,comparison=pp$contrast,estimate=pp$estimate,SE=pp$SE,df=pp$df,lower=pp$lower.CL,upper=pp$upper.CL,p_value=pp$p.value,p_bh=p.adjust(pp$p.value,'BH',n=3),n_donors=nlevels(z$case_id),n_objects=nrow(z))
  aa<-as.data.frame(anova(fit,ddf='Satterthwaite'));term<-'exposure_z:donor_status'
  gr<-data.frame(model=name,source=label,role=role,comparison='ND vs Aab+ vs T1D: adjusted proximity slope',test='Satterthwaite F test of disease-by-proximity interaction',p_value=aa[term,'Pr(>F)'],statistic=aa[term,'F value'],num_df=aa[term,'NumDF'],den_df=aa[term,'DenDF'],n_donors=nlevels(z$case_id),n_objects=nrow(z),status='estimated')
  list(slopes=sr,pairs=pr,global=gr)
 }
 if(!is.null(chosen)){
  extracted<-extract_fit(chosen,paste0('random_slope_',chosen_covariance),'declared_model')
  slopes[[name]]<-extracted$slopes;pairs_out[[name]]<-extracted$pairs;global[[name]]<-extracted$global
  cf<-as.data.frame(coef(summary(chosen)));cf$term<-rownames(cf);cf$model<-name;fixed[[name]]<-cf
 }else{
  dz<-dc[dc$estimable,];dz$donor_status<-factor(dz$donor_status,levels=c('ND','Aab+','T1D'))
  replicated<-all(table(dz$donor_status)>=2)
  kw<-if(replicated)kruskal.test(estimate~donor_status,dz)else NULL
  global[[name]]<-data.frame(model=name,source='independent_donor_regression',role='declared_fallback',comparison='ND vs Aab+ vs T1D: adjusted proximity slope',test='Kruskal-Wallis on independent donor coefficients',p_value=if(is.null(kw))NA else kw$p.value,statistic=if(is.null(kw))NA else unname(kw$statistic),num_df=2,den_df=NA,n_donors=nrow(dz),n_objects=sum(dz$objects),status=if(replicated)'estimated'else'unestimable: fewer than two donor coefficients in a state')
  sl<-list();pa<-list()
  for(g in levels(d$donor_status)){
   v<-dz$estimate[dz$donor_status==g];tt<-if(length(v)>=2&&sd(v)>1e-12)t.test(v)else NULL
   sl[[g]]<-data.frame(model=name,source='independent_donor_regression',role='declared_fallback',donor_status=g,estimate=if(length(v))mean(v)else NA,SE=if(length(v)>=2)sd(v)/sqrt(length(v))else NA,df=length(v)-1,lower=if(is.null(tt))NA else tt$conf.int[1],upper=if(is.null(tt))NA else tt$conf.int[2],p_value=if(is.null(tt))NA else tt$p.value,n_donors=length(v),n_objects=sum(dz$objects[dz$donor_status==g]))
  }
  ss<-bind(sl);ss$p_bh<-p.adjust(ss$p_value,'BH',n=3);slopes[[name]]<-ss
  for(groups in combn(levels(d$donor_status),2,simplify=FALSE)){
   a<-dz$estimate[dz$donor_status==groups[1]];b<-dz$estimate[dz$donor_status==groups[2]]
   p<-if(length(a)>=2&&length(b)>=2)suppressWarnings(wilcox.test(a,b,exact=!anyDuplicated(c(a,b)),correct=TRUE)$p.value)else NA
   pa[[paste(groups,collapse='_')]]<-data.frame(model=name,source='independent_donor_regression',role='declared_fallback',comparison=paste(groups,collapse=' - '),estimate=if(length(a)&&length(b))mean(a)-mean(b)else NA,SE=NA,df=NA,lower=NA,upper=NA,p_value=p,n_donors=length(a)+length(b),n_objects=sum(dz$objects[dz$donor_status%in%groups]))
  }
  pp<-bind(pa);pp$p_bh<-p.adjust(pp$p_value,'BH',n=3);pairs_out[[name]]<-pp
  # This familiar model is separately labelled, never promoted by its p value.
  ri<-lmer(as.formula(paste(ff,'+(1|case_id)')),z,REML=TRUE,control=lmerControl(optimizer='bobyqa'))
  saveRDS(ri,file.path(models,paste0('h7_',name,'_supporting_random_intercept.rds')))
  diag[[paste0(name,'_supporting_random_intercept')]]<-list(singular=isSingular(ri),optimizer_code=ri@optinfo$conv$opt,convergence_messages=ri@optinfo$conv$lme4$messages,variance=as.data.frame(VarCorr(ri)))
  ex<-extract_fit(ri,'random_intercept_only','supporting_sensitivity')
  slopes[[paste0(name,'_ri')]]<-ex$slopes;pairs_out[[paste0(name,'_ri')]]<-ex$pairs;global[[paste0(name,'_ri')]]<-ex$global
 }
 cat(name,'objects',nrow(z),'donors',nlevels(z$case_id),'zeroBeta',sum(z$n_beta_core==0),'selected',if(is.null(chosen))'donor_fallback'else chosen_covariance,'\n')
}
for(nm in c('slopes','pairs_out','global','donor_rows','eligibility','fixed','scaling'))write.csv(bind(get(nm)),file.path(models,paste0('h7_',nm,'.csv')),row.names=FALSE)
gg<-bind(global);gg$sensitivity_global_q<-NA_real_;ix<-gg$model!='primary_20um'&gg$role!='supporting_sensitivity';gg$sensitivity_global_q[ix]<-p.adjust(gg$p_value[ix],'BH',n=4);write.csv(gg,file.path(models,'h7_global.csv'),row.names=FALSE)
write_json(diag,file.path(models,'h7_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=16)
capture.output(sessionInfo(),file=file.path(models,'R_session.txt'))
