#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(lmerTest);library(emmeans);library(jsonlite)})
emm_options(lmerTest.limit=10000)
out<-'output/h5_h10_original_plan_20261026/h6'
raw<-read.csv(file.path(out,'h6_object_features.csv'),stringsAsFactors=FALSE)
raw$case_id<-factor(raw$case_id);raw$donor_status<-factor(raw$donor_status,levels=c('ND','Aab+','T1D'))
raw<-raw[tolower(as.character(raw$h6_primary_eligible))=='true',]
h7<-read.csv('output/h5_h10_original_plan_20261026/h7/object_neurovascular_features.csv')
raw<-merge(raw,h7[,c('islet_id','nerve_area_fraction','cd31_area_fraction','nerve_near_cd31_20um')],by='islet_id',all.x=TRUE,sort=FALSE)
base_controls<-c('beta_fraction','delta_fraction','other_endocrine_fraction','log_core_area','log_nearby_immune','log_core_immune')
spec<-data.frame(model=c('primary','alpha','delta','beta_density','contact_by_immune','raw_INS','corrected_INS','strong_hormone_margin','nerve_vessel_adjustment','boundary_contacts'),
 contact=c('combined_excess','alpha_excess','delta_excess',rep('combined_excess',6),'boundary_combined_excess'),
 outcome=c('ins_mean_beta_core','ins_mean_beta_core','ins_mean_beta_core','beta_density_mm2',
           'ins_mean_beta_core','ins_mean_beta_raw','ins_mean_beta_corrected','ins_mean_beta_margin4','ins_mean_beta_core','ins_mean_beta_core'))
spec<-rbind(spec,data.frame(model='nerve_vessel_proximity_adjustment',contact='combined_excess',outcome='ins_mean_beta_core'))
write.csv(spec,file.path(out,'model_specifications_before_fit.csv'),row.names=FALSE)
if(file.exists(file.path(out,'boundary_contact_measurements.csv'))){
 b<-read.csv(file.path(out,'boundary_contact_measurements.csv'))
 b$boundary_combined_excess<-ifelse(tolower(as.character(b$combined_spatially_informative))=='true',b$combined_excess,NA_real_)
 raw<-merge(raw,b[,c('islet_id','boundary_combined_excess')],by='islet_id',all.x=TRUE,sort=FALSE)
}else raw$boundary_combined_excess<-NA_real_
slopes<-list();contrasts<-list();tests<-list();coefs<-list();diags<-list();donors<-list();donortests<-list();donorpairs<-list();donor_direction<-list()
for(k in seq_len(nrow(spec))){
 name<-spec$model[k];contact<-spec$contact[k];outcome<-spec$outcome[k]
 controls<-base_controls
 if(name=='nerve_vessel_adjustment')controls<-c(controls,'nerve_area_fraction','cd31_area_fraction')
 if(name=='nerve_vessel_proximity_adjustment')controls<-c(controls,'nerve_area_fraction','cd31_area_fraction','nerve_near_cd31_20um')
 d<-droplevels(raw[is.finite(raw[[contact]])&is.finite(raw[[outcome]])&apply(raw[,controls],1,function(x)all(is.finite(x))),])
 if(name=='alpha')d<-droplevels(d[d$n_alpha_nodes>=2&tolower(as.character(d$alpha_spatially_informative))=='true',])
 if(name=='delta')d<-droplevels(d[d$n_delta_nodes>=2&tolower(as.character(d$delta_spatially_informative))=='true',])
 if(nrow(d)<10||nlevels(d$case_id)<3){cat('Not estimable:',name,'objects=',nrow(d),'\n');next}
 d$y_z<-as.numeric(scale(d[[outcome]]))
 zvars<-c('contact'=contact,setNames(controls,controls))
 for(v in names(zvars))d[[paste0(v,'_z')]]<-ave(d[[zvars[[v]]]],d$case_id,FUN=function(x)if(sd(x)<1e-12)rep(0,length(x))else as.numeric(scale(x)))
 zcontrols<-paste0(controls,'_z')
 write.csv(d,file.path(out,paste0(name,'_analysis_objects.csv')),row.names=FALSE)
 # Fixed full covariate set: donors with insufficient information remain in mixed
 # models but are explicitly unavailable in the independently estimated slopes.
 dr<-list();dfits<-list()
 for(case in levels(d$case_id)){
   dd<-droplevels(d[d$case_id==case,]);cc<-zcontrols[vapply(zcontrols,function(v)sd(dd[[v]])>1e-12,logical(1))]
   form<-as.formula(paste('y_z ~ contact_z +',paste(cc,collapse=' + ')))
   if(name=='contact_by_immune')form<-update(form,.~.+contact_z:log_nearby_immune_z)
   fit<-lm(form,data=dd);term<-if(name=='contact_by_immune')'contact_z:log_nearby_immune_z'else'contact_z'
   valid<-fit$rank==ncol(model.matrix(fit))&&df.residual(fit)>0&&is.finite(coef(fit)[term])
   rr<-data.frame(model=name,case_id=case,donor_status=as.character(dd$donor_status[1]),objects=nrow(dd),
     parameters=ncol(model.matrix(fit)),rank=fit$rank,residual_df=df.residual(fit),estimable=valid,
     standardized_beta=NA_real_,SE=NA_real_,lower=NA_real_,upper=NA_real_,formula=paste(deparse(form),collapse=' '))
   if(valid){r<-coef(summary(fit))[term,];ci<-confint(fit,term);rr$standardized_beta<-r[1];rr$SE<-r[2];rr$lower<-ci[1];rr$upper<-ci[2]}
   dr[[case]]<-rr;dfits[[case]]<-fit
 }
 dx<-do.call(rbind,dr);donors[[name]]<-dx;saveRDS(dfits,file.path(out,paste0(name,'_donor_fits.rds')))
 dx<-dx[dx$estimable,]
 if(nrow(dx)>=3){
   groups_supported<-all(table(dx$donor_status)>=2)
   kw<-if(groups_supported)kruskal.test(standardized_beta~donor_status,data=dx)else NULL
   donortests[[name]]<-data.frame(model=name,test='Kruskal-Wallis',p_value=if(is.null(kw))NA_real_ else kw$p.value,
      statistic=if(is.null(kw))NA_real_ else unname(kw$statistic),donors=nrow(dx),
      reason=if(groups_supported)'All states have replicated donor coefficients'else'T1D has only one independently estimable donor coefficient')
   pp<-list()
   for(pair in combn(c('ND','Aab+','T1D'),2,simplify=FALSE)){
     a<-dx$standardized_beta[dx$donor_status==pair[1]];b<-dx$standardized_beta[dx$donor_status==pair[2]]
     if(length(a)>=2&&length(b)>=2){w<-wilcox.test(a,b,exact=!anyDuplicated(c(a,b)),correct=TRUE);p<-w$p.value}else p<-NA_real_
     pp[[paste(pair,collapse='_')]]<-data.frame(model=name,group_1=pair[1],group_2=pair[2],n1=length(a),n2=length(b),p_value=p)
   }
   pp<-do.call(rbind,pp);pp$p_bh<-p.adjust(pp$p_value,'BH',n=3);donorpairs[[name]]<-pp
   for(group in c('All donors','ND','Aab+','T1D')){
     v<-if(group=='All donors')dx$standardized_beta else dx$standardized_beta[dx$donor_status==group]
     if(!length(v))next
     w<-wilcox.test(v,mu=0,exact=!anyDuplicated(abs(v))&&!any(v==0),correct=TRUE)
     tt<-if(length(v)>1)t.test(v,mu=0)else NULL
     donor_direction[[paste(name,group)]]<-data.frame(model=name,group=group,donors=length(v),mean_beta=mean(v),median_beta=median(v),
       mean_lower=if(is.null(tt))NA_real_ else tt$conf.int[1],mean_upper=if(is.null(tt))NA_real_ else tt$conf.int[2],
       positive_donors=sum(v>0),negative_donors=sum(v<0),signed_rank_p=if(length(v)<2)NA_real_ else w$p.value)
   }
 }
 for(association in c('common','disease_interaction')){
   fixed<-paste('y_z ~ contact_z',if(association=='common')'+ donor_status'else'* donor_status','+',paste(zcontrols,collapse=' + '))
   if(name=='contact_by_immune')fixed<-paste(fixed,'+ contact_z:log_nearby_immune_z')
   for(covariance in c('correlated','independent','intercept_only')){
     if(covariance!='correlated'&&!singular)next
     random<-switch(covariance,correlated='(1+contact_z|case_id)',independent='(1|case_id)+(0+contact_z|case_id)',intercept_only='(1|case_id)')
     key<-paste(name,association,covariance,sep='_');form<-as.formula(paste(fixed,'+',random))
     warnings<-character()
     fit<-withCallingHandlers(lmer(form,data=d,REML=TRUE,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
     singular<-isSingular(fit);saveRDS(fit,file.path(out,paste0(key,'.rds')))
     aa<-as.data.frame(anova(fit,type=3,ddf='Satterthwaite'));aa$term<-rownames(aa);aa$model<-key;tests[[key]]<-aa
     cc<-as.data.frame(coef(summary(fit)));cc$term<-rownames(cc);cc$model<-key;cc$lower<-cc$Estimate-qt(.975,cc$df)*cc$`Std. Error`;cc$upper<-cc$Estimate+qt(.975,cc$df)*cc$`Std. Error`;coefs[[key]]<-cc
     if(association=='disease_interaction'){
       tr<-emtrends(fit,~donor_status,var='contact_z',lmer.df='satterthwaite')
       ss<-as.data.frame(summary(tr,infer=c(TRUE,TRUE),adjust='none'));ss$model<-key;slopes[[key]]<-ss
       pp<-as.data.frame(summary(pairs(tr,adjust='none'),infer=c(TRUE,TRUE)));pp$p_bh<-p.adjust(pp$p.value,'BH');pp$model<-key;contrasts[[key]]<-pp
     }
     diags[[key]]<-list(objects=nrow(d),donors=nlevels(d$case_id),singular=singular,optimizer_code=fit@optinfo$conv$opt,
       final_messages=fit@optinfo$conv$lme4$messages,warnings=unique(warnings),random_covariance=as.data.frame(VarCorr(fit)),
       min_hessian_eigenvalue=min(eigen(fit@optinfo$derivs$Hessian,symmetric=TRUE,only.values=TRUE)$values),formula=paste(deparse(form),collapse=' '))
     cat(key,'singular=',singular,'n=',nrow(d),'\n');print(cc[cc$term%in%c('contact_z','contact_z:log_nearby_immune_z'),]);if(association=='disease_interaction')print(aa[aa$term=='contact_z:donor_status',])
   }
 }
}
write.csv(do.call(rbind,donors),file.path(out,'donor_regression_eligibility.csv'),row.names=FALSE)
write.csv(do.call(rbind,donortests),file.path(out,'donor_disease_global_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,donorpairs),file.path(out,'donor_disease_pairwise_tests.csv'),row.names=FALSE)
dd<-do.call(rbind,donor_direction);dd$within_state_p_bh<-NA_real_
for(name in unique(dd$model)){
 ix<-dd$model==name&dd$group!='All donors'
 dd$within_state_p_bh[ix]<-p.adjust(dd$signed_rank_p[ix],'BH',n=3)
}
write.csv(dd,file.path(out,'donor_direction_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,slopes),file.path(out,'disease_slopes.csv'),row.names=FALSE)
write.csv(do.call(rbind,contrasts),file.path(out,'disease_slope_comparisons.csv'),row.names=FALSE)
write.csv(do.call(rbind,tests),file.path(out,'model_anova.csv'),row.names=FALSE)
write.csv(do.call(rbind,coefs),file.path(out,'model_coefficients.csv'),row.names=FALSE)
write_json(diags,file.path(out,'model_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
audit<-lapply(names(diags),function(key){
 fit<-readRDS(file.path(out,paste0(key,'.rds')))
 data.frame(model=key,objects=nobs(fit),design_rank=qr(getME(fit,'X'))$rank,
   design_columns=ncol(getME(fit,'X')),singular=isSingular(fit),optimizer_code=fit@optinfo$conv$opt,
   min_hessian_eigenvalue=min(eigen(fit@optinfo$derivs$Hessian,symmetric=TRUE,only.values=TRUE)$values))
})
write.csv(do.call(rbind,audit),file.path(out,'final_rds_convergence_audit.csv'),row.names=FALSE)
