#!/usr/bin/env Rscript
# Direct disease-state comparisons of H2. Existing pixels/features are unchanged.
# Local precedent: Barlow2025 random donor intercept/slope with standardized
# measurements; Bodenmiller2025 density mixed models and multiple-comparison FDR.
options(warn=1)
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(emmeans))
out <- 'output/h1_h4_disease_comparisons_20261026/h2'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
input <- 'output/h1_h4_literature_methods_20261025/h2/h2_analysis_objects.csv'
d <- read.csv(input, stringsAsFactors=FALSE)
d$case_id <- factor(d$case_id)
d$donor_status <- factor(d$donor_status, levels=c('ND','Aab+','T1D'))
d$object_class <- factor(d$object_class)
d$core_density <- d$n_immune_core * 1e6 / d$core_area_um2
d$log_density <- log1p(d$core_density)
d$peri_density_log <- log1p(d$n_immune_peri * 1e6 / d$peri_area_um2)
d$log_core_area <- log(d$core_area_um2)
d$log_perimeter <- log(d$core_perimeter_um)
d$log_edge_distance <- log1p(pmax(d$core_distance_to_tissue_edge_um,0))
controls <- c('peri_density_log','log_core_area','log_perimeter','log_edge_distance','capillary_core_area_fraction')
within_z <- function(x,g) ave(x,g,FUN=function(v) if(sd(v)>1e-12) (v-mean(v))/sd(v) else rep(0,length(v)))
donor_scaling <- list()
for (v in controls) {
  donor_scaling[[v]] <- data.frame(case_id=levels(d$case_id),variable=v,
                                   mean=tapply(d[[v]],d$case_id,mean),sd=tapply(d[[v]],d$case_id,sd))
  d[[v]] <- within_z(d[[v]],d$case_id)
}
control_formula <- paste(c(controls,'object_class'), collapse=' + ')
stopifnot(nrow(d)==2937, nlevels(d$case_id)==20, !anyDuplicated(d$islet_id))
models <- list(); diagnostics <- list(); slopes <- list(); comparisons <- list()
means <- list(); mean_comparisons <- list(); tests <- list(); scaling <- list()
for (width in c(2,5)) for (outcome in c('core_density','log_density')) {
  name <- paste0('band',width,'_',outcome)
  z <- d
  exposure <- z[[paste0('coliv_coverage_',width,'um')]]
  keep <- is.finite(exposure) & is.finite(z[[outcome]])
  z <- droplevels(z[keep,]); exposure <- exposure[keep]
  ymean <- mean(z[[outcome]]); ysd <- sd(z[[outcome]])
  xmean <- mean(exposure); xsd <- sd(exposure)
  z$coverage_z <- within_z(exposure,z$case_id)
  donor_scaling[[name]] <- data.frame(case_id=levels(z$case_id),variable=paste0('coverage_',width,'um'),
                                      mean=tapply(exposure,z$case_id,mean),sd=tapply(exposure,z$case_id,sd))
  z$response_z <- (z[[outcome]]-ymean)/ysd
  scaling[[name]] <- data.frame(model=name,width_um=width,outcome=outcome,
                                coverage_mean=xmean,coverage_sd=xsd,response_mean=ymean,response_sd=ysd)
  fixed <- paste('response_z ~ coverage_z * donor_status +',control_formula)
  attempts <- list(); chosen <- NULL
  for (structure in c('correlated_donor_intercept_slope','independent_donor_intercept_slope')) {
    random <- if(structure=='correlated_donor_intercept_slope') '(1 + coverage_z | case_id)' else '(1 + coverage_z || case_id)'
    form <- as.formula(paste(fixed,'+',random))
    caught <- character()
    fit <- withCallingHandlers(lmer(form,z,REML=TRUE,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),
                               warning=function(w){caught <<- c(caught,conditionMessage(w));invokeRestart('muffleWarning')})
    variance <- as.data.frame(VarCorr(fit))
    slope_sd <- variance$sdcor[!is.na(variance$var1) & variance$var1=='coverage_z' & is.na(variance$var2)]
    convergence <- fit@optinfo$conv$lme4$messages
    attempts[[structure]] <- list(singular=isSingular(fit),slope_sd=slope_sd,warnings=c(caught,convergence),variance=variance)
    saveRDS(fit,file.path(out,paste0(name,'_',structure,'.rds')))
    if(!isSingular(fit) && length(slope_sd)==1 && slope_sd>1e-6 && length(convergence)==0) {
      chosen <- fit; chosen_structure <- structure; break
    }
  }
  diagnostics[[name]] <- attempts
  if(is.null(chosen)) {
    saveRDS(diagnostics,file.path(out,'model_diagnostics.rds'))
    capture.output(str(diagnostics),file=file.path(out,'model_diagnostics.txt'))
    message(paste(name,'has no nonsingular donor random-slope model; direct disease inference has not been substituted with object-level inference'))
    next
  }
  models[[name]] <- chosen
  a <- as.data.frame(anova(chosen)); a$term <- rownames(a)
  a$model <- name; a$width_um <- width; a$outcome <- outcome
  tests[[name]] <- a
  tr <- emtrends(chosen,~donor_status,var='coverage_z',lmer.df='satterthwaite')
  st <- as.data.frame(summary(tr,infer=c(TRUE,TRUE),adjust='none'))
  st$model <- name; st$width_um <- width; st$outcome <- outcome; st$random_structure <- chosen_structure
  slopes[[name]] <- st
  cp <- as.data.frame(summary(pairs(tr),infer=c(TRUE,TRUE),adjust='none'))
  cp$p_BH <- p.adjust(cp$p.value,method='BH')
  cp$model <- name; cp$width_um <- width; cp$outcome <- outcome
  comparisons[[name]] <- cp
  em <- emmeans(chosen,~donor_status,at=list(coverage_z=0),weights='proportional',lmer.df='satterthwaite')
  mt <- as.data.frame(summary(em,infer=c(TRUE,TRUE),adjust='none'))
  mt$adjusted_mean_original_units <- ymean+ysd*mt$emmean
  mt$lower_original_units <- ymean+ysd*mt$lower.CL
  mt$upper_original_units <- ymean+ysd*mt$upper.CL
  mt$model <- name; mt$width_um <- width; mt$outcome <- outcome
  means[[name]] <- mt
  mc <- as.data.frame(summary(pairs(em),infer=c(TRUE,TRUE),adjust='none'))
  mc$p_BH <- p.adjust(mc$p.value,method='BH')
  mc$difference_original_units <- ysd*mc$estimate
  mc$lower_original_units <- ysd*mc$lower.CL
  mc$upper_original_units <- ysd*mc$upper.CL
  mc$model <- name; mc$width_um <- width; mc$outcome <- outcome
  mean_comparisons[[name]] <- mc
  print(name); print(a[a$term=='coverage_z:donor_status',]); print(cp)
}
write.csv(do.call(rbind,tests),file.path(out,'model_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,slopes),file.path(out,'standardized_slopes.csv'),row.names=FALSE)
write.csv(do.call(rbind,comparisons),file.path(out,'disease_slope_comparisons.csv'),row.names=FALSE)
write.csv(do.call(rbind,means),file.path(out,'disease_adjusted_means.csv'),row.names=FALSE)
write.csv(do.call(rbind,mean_comparisons),file.path(out,'disease_mean_comparisons.csv'),row.names=FALSE)
write.csv(do.call(rbind,scaling),file.path(out,'standardization.csv'),row.names=FALSE)
write.csv(do.call(rbind,donor_scaling),file.path(out,'donor_standardization.csv'),row.names=FALSE)
saveRDS(models,file.path(out,'models.rds'))
saveRDS(diagnostics,file.path(out,'model_diagnostics.rds'))
capture.output(str(diagnostics),file=file.path(out,'model_diagnostics.txt'))

# The primary raw-density random-slope covariance is singular. Preserve those
# fits above, and compare one adjusted coefficient per donor with the published
# donor-level rank-test approach. This rule depends on model identifiability,
# not on the resulting group-comparison p-values.
if (!('band2_core_density' %in% names(models))) {
  d$response_z <- as.numeric(scale(d$core_density))
  d$coverage_z <- within_z(d$coliv_coverage_2um,d$case_id)
  donor_fits <- list(); donor_rows <- list()
  for (case in levels(d$case_id)) {
    z <- droplevels(d[d$case_id==case,])
    varying <- controls[sapply(z[controls],function(v) sd(v)>1e-12)]
    terms <- c('coverage_z',varying,if(nlevels(z$object_class)>1)'object_class')
    form <- reformulate(terms,response='response_z')
    m <- lm(form,z)
    if(m$rank<ncol(model.matrix(m)) || !is.finite(coef(m)[['coverage_z']]))
      stop(paste('Donor coefficient not identifiable:',case))
    cf <- coef(summary(m))['coverage_z',]
    ci <- confint(m,'coverage_z')
    donor_fits[[case]] <- m
    donor_rows[[case]] <- data.frame(case_id=case,donor_status=as.character(z$donor_status[1]),
      objects=nrow(z),core_immune_cells=sum(z$n_immune_core),standardized_beta=cf[['Estimate']],
      within_donor_se=cf[['Std. Error']],within_donor_lower=ci[1],within_donor_upper=ci[2],
      residual_df=df.residual(m),design_rank=m$rank,design_columns=ncol(model.matrix(m)))
  }
  donor <- do.call(rbind,donor_rows)
  donor$donor_status <- factor(donor$donor_status,levels=c('ND','Aab+','T1D'))
  kw <- kruskal.test(standardized_beta~donor_status,data=donor)
  overall <- data.frame(method=kw$method,statistic=unname(kw$statistic),df=unname(kw$parameter),
                        p=kw$p.value,donors=nrow(donor),endpoint='2um_raw_core_density_standardized_adjusted_donor_coefficient')
  pair_rows <- list()
  for (pair in combn(levels(donor$donor_status),2,simplify=FALSE)) {
    a <- donor$standardized_beta[donor$donor_status==pair[1]]
    b <- donor$standardized_beta[donor$donor_status==pair[2]]
    wt <- wilcox.test(a,b,exact=TRUE,alternative='two.sided')
    pair_rows[[paste(pair,collapse='_')]] <- data.frame(group1=pair[1],group2=pair[2],n1=length(a),n2=length(b),
      median1=median(a),median2=median(b),median_difference=median(a)-median(b),
      statistic=unname(wt$statistic),p=wt$p.value,method=wt$method)
  }
  pair_table <- do.call(rbind,pair_rows)
  pair_table$p_BH <- p.adjust(pair_table$p,method='BH')
  group_table <- do.call(rbind,lapply(levels(donor$donor_status),function(g) {
    v <- donor$standardized_beta[donor$donor_status==g]
    data.frame(donor_status=g,donors=length(v),median=median(v),q25=unname(quantile(v,.25)),q75=unname(quantile(v,.75)))
  }))
  write.csv(donor,file.path(out,'primary_donor_coefficients.csv'),row.names=FALSE)
  write.csv(overall,file.path(out,'primary_disease_overall_test.csv'),row.names=FALSE)
  write.csv(pair_table,file.path(out,'primary_disease_pairwise_tests.csv'),row.names=FALSE)
  write.csv(group_table,file.path(out,'primary_disease_summary.csv'),row.names=FALSE)
  saveRDS(donor_fits,file.path(out,'primary_donor_regressions.rds'))
  print('Primary donor-level disease comparisons'); print(overall); print(pair_table); print(group_table)
}
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
