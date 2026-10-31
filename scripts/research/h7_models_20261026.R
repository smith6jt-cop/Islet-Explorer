#!/usr/bin/env Rscript
# Approved observed-proximity H7 substitution; no spatial-null endpoint.
options(warn=1)
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(jsonlite))
out <- 'output/h5_h10_original_plan_20261026/h7'
d <- read.csv(file.path(out,'object_neurovascular_features.csv'),stringsAsFactors=FALSE)
d$case_id <- factor(d$case_id)
d$donor_status <- factor(d$donor_status,levels=c('ND','Aab+','T1D'))
d$neighborhood_complete_20 <- tolower(as.character(d$neighborhood_complete_20))=='true'
d$log_core_area <- log(d$core_area_um2)
d$log_perimeter <- log(d$core_perimeter_um)
d$log_nearby_immune_density <- log1p(1e6*d$n_immune_nearby/(d$core_area_um2+d$peri_area_um2))
controls <- c('nerve_area_fraction','cd31_area_fraction','log_core_area','log_perimeter','log_nearby_immune_density')
spec <- data.frame(model=c('primary_20um','distance_10um','distance_5um','overlap','GAP43_20um','raw_INS_20um','corrected_INS_20um'),
  exposure=c('nerve_near_cd31_20um','nerve_near_cd31_10um','nerve_near_cd31_5um','nerve_vessel_overlap','gap43_nerve_near_cd31_20um','nerve_near_cd31_20um','nerve_near_cd31_20um'),
  outcome=c(rep('ins_mean_beta_core',5),'ins_mean_beta_raw','ins_mean_beta_corrected'),
  role=c('primary',rep('sensitivity',6)))
write.csv(spec,file.path(out,'inference_model_specifications.csv'),row.names=FALSE)
write_json(list(approval='User approved adjusted observed-proximity regression replacing excess-over-spatial-null alignment.',
  endpoint='Within-Beta mean INS; global standardization, no outcome log transformation',
  exposure='Within-donor standardized observed proximity; primary 20 micrometres',
  covariates=controls,covariate_scaling='Within-donor z on the complete-case set for each model; no required covariate omitted',
  fixed='response_z ~ exposure_z * donor_status + nerve_area_fraction_z + cd31_area_fraction_z + log_core_area_z + log_perimeter_z + log_nearby_immune_density_z',
  donor_replication='Random donor intercept and exposure slope; Satterthwaite inference from lmerTest',
  covariance_rule='Correlated intercept/slope first; independent components only when correlated covariance is singular or not converged. Never choose by p value.',
  fallback='If neither fit is nonsingular/converged, retain diagnostic models and calculate independent full-covariate donor regression coefficients only if exposure varies, design is full rank and residual degrees of freedom > 0; disease tests on estimable donor coefficients only.',
  hypothesis_tests='Conditional disease-specific slopes versus zero and all three disease slope differences. Separate BH families of three tests within each specification. Equal-weight mean of three disease slopes additionally describes the common conditional direction.',
  meaningful_effect='No minimum important effect was specified. Original-plan branch judgment is not assigned.',
  local_precedent=c('literature/Barlow2025.pdf: Methods random donor intercept/slope, within-donor predictor standardization, global outcome standardization, Satterthwaite inference',
  'literature/Bodenmiller2025.pdf: donor-level Wilcoxon comparisons and FDR; literature/Richardson2023.pdf: Kruskal-Wallis',
  'literature/Heide.pdf: one-sample t test and confidence intervals; donor-slope mean is a conditional fallback, not a correlation')),
  file.path(out,'inference_analysis_plan.json'),pretty=TRUE,auto_unbox=TRUE)
within_z <- function(x,g) ave(x,g,FUN=function(v) if(sd(v)>1e-12)(v-mean(v))/sd(v) else rep(0,length(v)))
bind <- function(x) if(length(x))do.call(rbind,x) else data.frame()
fits <- list(); diagnostics <- list(); slopes <- list(); pairs_out <- list(); common <- list(); tests <- list(); fixed <- list(); eligibility <- list(); scaling <- list(); donor_rows <- list(); fallback_tests <- list(); fallback_pairs <- list(); fallback_means <- list(); design <- list()
for(i in seq_len(nrow(spec))) {
  name <- spec$model[i]; exposure <- spec$exposure[i]; outcome <- spec$outcome[i]
  use <- is.finite(d[[exposure]]) & is.finite(d[[outcome]]) & d$neighborhood_complete_20
  for(v in controls)use <- use & is.finite(d[[v]])
  z <- droplevels(d[use,]); z$exposure_z <- within_z(z[[exposure]],z$case_id)
  ym <- mean(z[[outcome]]); ys <- sd(z[[outcome]])
  z$response_z <- (z[[outcome]]-ym)/ys
  for(v in controls)z[[paste0(v,'_z')]] <- within_z(z[[v]],z$case_id)
  write.csv(z,file.path(out,paste0(name,'_model_objects.csv')),row.names=FALSE)
  eligibility[[name]] <- bind(lapply(levels(d$case_id),function(ca){
    all <- d[d$case_id==ca,]; zz <- z[z$case_id==ca,]
    data.frame(model=name,case_id=ca,donor_status=as.character(all$donor_status[1]),h5_objects=nrow(all),beta_objects=sum(is.finite(all[[outcome]])),
      matched_objects=nrow(zz),exposure_sd=if(nrow(zz)>1)sd(zz[[exposure]]) else NA_real_,
      exposure_varies=nrow(zz)>1 && sd(zz[[exposure]])>1e-12,
      absent_nerve_or_corroborating_signal=sum(is.finite(all[[outcome]]) & all$neighborhood_complete_20 & !is.finite(all[[exposure]])),
      incomplete_20um=sum(is.finite(all[[outcome]]) & !all$neighborhood_complete_20))}))
  scaling[[name]] <- data.frame(model=name,outcome=outcome,mean=ym,sd=ys,objects=nrow(z),donors=nlevels(z$case_id))
  for(ca in levels(z$case_id))for(v in c(exposure,controls)){
    xx <- z[z$case_id==ca,v]
    design[[paste(name,ca,v)]] <- data.frame(model=name,case_id=ca,variable=v,mean=mean(xx),sd=sd(xx))
  }
  control_names <- paste0(controls,'_z')
  ff <- paste('response_z ~ exposure_z * donor_status +',paste(control_names,collapse=' + '))
  chosen <- NULL; attempts <- list()
  for(structure in c('correlated','independent')){
    random <- if(structure=='correlated')'(1 + exposure_z | case_id)' else '(1 + exposure_z || case_id)'
    warned <- character(); form <- as.formula(paste(ff,'+',random))
    fit <- withCallingHandlers(lmer(form,z,REML=TRUE,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),warning=function(w){warned<<-c(warned,conditionMessage(w));invokeRestart('muffleWarning')})
    conv <- fit@optinfo$conv$lme4$messages
    vc <- as.data.frame(VarCorr(fit)); singular <- isSingular(fit)
    attempts[[structure]] <- list(formula=paste(deparse(form),collapse=' '),singular=singular,convergence_messages=conv,warnings=warned,variance=vc)
    saveRDS(fit,file.path(out,paste0(name,'_',structure,'.rds')))
    if(!singular && length(conv)==0){chosen <- fit; chosen_structure <- structure; break}
  }
  diagnostics[[name]] <- attempts
  # Fully adjusted individual donor coefficients are exported even where the
  # main mixed model is valid; inadequate donors remain explicit NA rows.
  donor_fit_list <- list(); donor_current <- list()
  for(ca in levels(z$case_id)){
    zz <- z[z$case_id==ca,]
    fm <- lm(reformulate(c('exposure_z',control_names),response='response_z'),zz)
    varying <- sd(zz$exposure_z)>1e-12
    estimable <- varying && fm$rank==ncol(model.matrix(fm)) && df.residual(fm)>0 && is.finite(coef(fm)['exposure_z'])
    cf <- if(estimable)coef(summary(fm))['exposure_z',] else rep(NA_real_,4)
    ci <- if(estimable)confint(fm,'exposure_z') else matrix(c(NA,NA),nrow=1)
    donor_current[[ca]] <- data.frame(model=name,case_id=ca,donor_status=as.character(zz$donor_status[1]),objects=nrow(zz),
      exposure_varies=varying,design_rank=fm$rank,design_columns=ncol(model.matrix(fm)),residual_df=df.residual(fm),
      estimable=estimable,standardized_beta=cf[1],SE=cf[2],lower=ci[1],upper=ci[2],
      reason=if(!varying)'constant exposure' else if(fm$rank<ncol(model.matrix(fm)))'insufficient independent observations for all required covariates' else if(df.residual(fm)<=0)'no residual degrees of freedom' else 'estimable')
    donor_fit_list[[ca]] <- fm
  }
  dc <- bind(donor_current); donor_rows[[name]] <- dc
  saveRDS(donor_fit_list,file.path(out,paste0(name,'individual_donor_models.rds')))
  if(is.null(chosen)){
    dz <- dc[dc$estimable,]; dz$donor_status <- factor(dz$donor_status,levels=c('ND','Aab+','T1D'))
    if(length(unique(dz$donor_status))>=2){
      kw <- kruskal.test(standardized_beta~donor_status,dz)
      fallback_tests[[name]] <- data.frame(model=name,method=kw$method,statistic=unname(kw$statistic),df=unname(kw$parameter),p=kw$p.value,donors=nrow(dz))
    }
    pr <- list(); gm <- list()
    for(pa in combn(levels(dz$donor_status),2,simplify=FALSE)){
      a <- dz$standardized_beta[dz$donor_status==pa[1]]; b <- dz$standardized_beta[dz$donor_status==pa[2]]
      p <- if(length(a)>=2 && length(b)>=2)wilcox.test(a,b,exact=TRUE)$p.value else NA_real_
      pr[[paste(pa,collapse='_')]] <- data.frame(model=name,group1=pa[1],group2=pa[2],n1=length(a),n2=length(b),median1=median(a),median2=median(b),p=p)
    }
    pt <- bind(pr); pt$p_BH <- p.adjust(pt$p,'BH'); fallback_pairs[[name]] <- pt
    for(g in levels(dz$donor_status)){
      a <- dz$standardized_beta[dz$donor_status==g]
      tm <- if(length(a)>=2)t.test(a,mu=0) else NULL
      gm[[g]] <- data.frame(model=name,donor_status=g,donors=length(a),mean=mean(a),median=median(a),
        lower=if(!is.null(tm))tm$conf.int[1] else NA,upper=if(!is.null(tm))tm$conf.int[2] else NA,p=if(!is.null(tm))tm$p.value else NA,
        note='Mean adjusted donor coefficient; one-sample t interval/test conditional on estimable donors')
    }
    gt <- bind(gm); gt$p_BH <- p.adjust(gt$p,'BH'); fallback_means[[name]] <- gt
    message(name,': both random-slope fits singular/nonconverged; fallback estimable donors ',paste(table(dz$donor_status),collapse='/'))
    next
  }
  fits[[name]] <- chosen
  tt <- as.data.frame(anova(chosen,ddf='Satterthwaite'));tt$term<-rownames(tt);tt$model<-name;tests[[name]]<-tt
  tr <- emtrends(chosen,~donor_status,var='exposure_z',lmer.df='satterthwaite')
  ss <- as.data.frame(summary(tr,infer=c(TRUE,TRUE),adjust='none'));ss$p_BH<-p.adjust(ss$p.value,'BH');ss$model<-name;ss$random_structure<-chosen_structure;slopes[[name]]<-ss
  pp <- as.data.frame(summary(pairs(tr),infer=c(TRUE,TRUE),adjust='none'));pp$p_BH<-p.adjust(pp$p.value,'BH');pp$model<-name;pairs_out[[name]]<-pp
  cc <- as.data.frame(summary(contrast(tr,list(equal_disease_weight_mean=c(1/3,1/3,1/3))),infer=c(TRUE,TRUE),adjust='none'));cc$model<-name;common[[name]]<-cc
  fx <- as.data.frame(coef(summary(chosen)));fx$term<-rownames(fx);fx$model<-name;fixed[[name]]<-fx
  cat('\n',name,' ',chosen_structure,'\n',sep='');print(ss);print(pp)
}
for(nm in c('slopes','pairs_out','common','tests','fixed','eligibility','scaling','donor_rows','fallback_tests','fallback_pairs','fallback_means','design')){
  filenames <- c(slopes='conditional_slopes.csv',pairs_out='disease_slope_comparisons.csv',common='equal_disease_mean_slopes.csv',tests='mixed_model_tests.csv',fixed='all_fixed_coefficients.csv',eligibility='model_donor_eligibility.csv',scaling='outcome_standardization.csv',donor_rows='individual_donor_coefficients.csv',fallback_tests='fallback_disease_tests.csv',fallback_pairs='fallback_disease_pairs.csv',fallback_means='fallback_conditional_means.csv',design='donor_predictor_standardization.csv')
  write.csv(bind(get(nm)),file.path(out,filenames[[nm]]),row.names=FALSE)
}
saveRDS(fits,file.path(out,'valid_inference_models.rds'));saveRDS(diagnostics,file.path(out,'inference_model_diagnostics.rds'))
write_json(diagnostics,file.path(out,'inference_model_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
capture.output(sessionInfo(),file=file.path(out,'inference_R_session.txt'))
