#!/usr/bin/env Rscript
# Direct disease-state comparison of within-donor standardized H1 associations.
# Local precedent: literature/Barlow2025.pdf p22, within-donor z-standardized
# predictors and donor random intercepts/slopes; Satterthwaite tests via lmerTest.
suppressPackageStartupMessages({library(lme4); library(lmerTest); library(emmeans)})
out <- 'output/h1_h4_disease_comparisons_20261026/h1'
source <- 'output/h1_h4_literature_methods_20261025/h1/analysis_objects.csv'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- read.csv(source, stringsAsFactors=FALSE)
stopifnot(nrow(d)==2937, length(unique(d$case_id))==20, !any(d$case_id==78))
d$case_id <- factor(d$case_id)
d$donor_status <- factor(d$donor_status, levels=c('ND','Aab+','T1D'))
d$object_class <- factor(d$object_class)
continuous <- c('immune_localization_fraction','log_area','log_area_ratio','log_nearby','log_field')
standardization <- list()
for (v in continuous) {
  output <- if (v=='immune_localization_fraction') 'exposure_z' else paste0(v,'_z')
  d[[output]] <- NA_real_
  for (case in levels(d$case_id)) {
    rows <- d$case_id==case
    mu <- mean(d[[v]][rows]); sdv <- sd(d[[v]][rows])
    stopifnot(is.finite(mu), is.finite(sdv))
    if (v=='immune_localization_fraction' && sdv==0) stop('No exposure variation in donor ',case)
    d[[output]][rows] <- if (sdv>0) (d[[v]][rows]-mu)/sdv else 0
    standardization[[length(standardization)+1]] <- data.frame(case_id=case,variable=v,
      donor_status=as.character(d$donor_status[which(rows)[1]]),objects=sum(rows),mean=mu,sd=sdv,
      zero_within_donor_variance=sdv==0)
  }
}
y_mean <- mean(d$beta_density_mm2); y_sd <- sd(d$beta_density_mm2)
d$y_z <- (d$beta_density_mm2-y_mean)/y_sd
write.csv(do.call(rbind,standardization),file.path(out,'within_donor_standardization.csv'),row.names=FALSE)
write.csv(data.frame(variable='beta_density_mm2',mean=y_mean,sd=y_sd,objects=nrow(d)),
          file.path(out,'global_outcome_standardization.csv'),row.names=FALSE)
write.csv(d,file.path(out,'analysis_objects_standardized.csv'),row.names=FALSE)

controls <- 'log_area_z + log_area_ratio_z + log_nearby_z + log_field_z + object_class'
fixed <- paste('y_z ~ exposure_z * donor_status +',controls)
diagnostics <- list(); groups <- list(); contrasts <- list(); interactions <- list()

fit_model <- function(name,random) {
  formula <- as.formula(paste(fixed,'+',random))
  warnings <- character()
  fit <- withCallingHandlers(lmerTest::lmer(formula,data=d,REML=TRUE,
    control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  saveRDS(fit,file.path(out,paste0(name,'.rds')))
  vc <- as.data.frame(VarCorr(fit))
  write.csv(vc,file.path(out,paste0(name,'_variance_components.csv')),row.names=FALSE)
  slope_var <- vc$vcov[which(vc$var1=='exposure_z' & is.na(vc$var2))]
  intercept_var <- vc$vcov[which(vc$var1=='(Intercept)' & is.na(vc$var2))]
  intercept_slope_correlation <- vc$sdcor[which(vc$var1=='(Intercept)' & vc$var2=='exposure_z')]
  if (!length(intercept_slope_correlation)) intercept_slope_correlation <- NA_real_
  hessian <- fit@optinfo$derivs$Hessian
  heig <- eigen(hessian,symmetric=TRUE,only.values=TRUE)$values
  gradient <- fit@optinfo$derivs$gradient
  scaled_gradient <- tryCatch(solve(chol(hessian),gradient),error=function(e)rep(NA_real_,length(gradient)))
  messages <- fit@optinfo$conv$lme4$messages
  diagnostics[[name]] <<- data.frame(model=name,objects=nobs(fit),donors=20,
    singular=isSingular(fit,tol=1e-4),optimizer_code=paste(fit@optinfo$conv$opt,collapse=';'),
    convergence_messages=paste(messages,collapse=' | '),warnings=paste(unique(warnings),collapse=' | '),
    minimum_hessian_eigenvalue=min(heig),maximum_absolute_gradient=max(abs(gradient)),
    maximum_absolute_scaled_gradient=max(abs(scaled_gradient)),
    random_intercept_variance=intercept_var,random_slope_variance=slope_var,
    intercept_slope_correlation=intercept_slope_correlation,
    residual_variance=sigma(fit)^2,formula=paste(deparse(formula),collapse=' '))
  a <- as.data.frame(anova(fit,type=3,ddf='Satterthwaite'))
  a$term <- rownames(a); a$model <- name
  write.csv(a,file.path(out,paste0(name,'_anova.csv')),row.names=FALSE)
  interaction <- a[a$term=='exposure_z:donor_status',]
  interactions[[name]] <<- interaction
  trend <- emtrends(fit,~donor_status,var='exposure_z',lmer.df='satterthwaite')
  g <- as.data.frame(summary(trend,infer=c(TRUE,TRUE),adjust='none'))
  g$model <- name
  g$donors <- as.integer(c(7,7,6)[match(g$donor_status,c('ND','Aab+','T1D'))])
  g$objects <- as.integer(c(792,1057,1088)[match(g$donor_status,c('ND','Aab+','T1D'))])
  groups[[name]] <<- g
  contrast <- as.data.frame(summary(pairs(trend,adjust='none'),infer=c(TRUE,TRUE),adjust='none'))
  contrast$p_bh <- p.adjust(contrast$p.value,method='BH')
  contrast$model <- name
  contrasts[[name]] <<- contrast
  print(diagnostics[[name]])
  print(interaction)
  print(g)
  print(contrast)
  fit
}

primary <- fit_model('correlated_donor_slopes','(1 + exposure_z | case_id)')
# Keep the requested correlated model as primary. If singular, show the prespecified
# independent intercept/slope sensitivity explicitly; never replace with intercept-only.
if (isSingular(primary,tol=1e-4) || length(primary@optinfo$conv$lme4$messages)>0) {
  sensitivity <- fit_model('independent_donor_slopes','(1 | case_id) + (0 + exposure_z | case_id)')
}

# The correlated covariance is singular, and independent random slopes reach zero.
# Use separately adjusted donor coefficients for the primary disease comparison.
# This fallback is triggered by covariance fit, never by significance.
donor_rows <- list(); donor_models <- list()
for (case in levels(d$case_id)) {
  z <- droplevels(d[d$case_id==case,])
  donor_controls <- c('log_area_z','log_area_ratio_z','log_nearby_z','log_field_z')
  donor_controls <- donor_controls[vapply(donor_controls,function(v)sd(z[[v]])>0,logical(1))]
  if (nlevels(z$object_class)>1) donor_controls <- c(donor_controls,'object_class')
  form <- as.formula(paste('y_z ~ exposure_z +',paste(donor_controls,collapse=' + ')))
  fit <- lm(form,data=z)
  coef <- coef(summary(fit))['exposure_z',]
  stopifnot(all(is.finite(coef)),!is.na(coef['Estimate']))
  ci <- confint(fit,'exposure_z',level=.95)
  donor_models[[case]] <- fit
  donor_rows[[case]] <- data.frame(case_id=case,donor_status=as.character(z$donor_status[1]),
    objects=nobs(fit),standardized_beta=coef['Estimate'],SE=coef['Std. Error'],
    lower=ci[1],upper=ci[2],residual_df=df.residual(fit),model_rank=fit$rank,
    full_rank=fit$rank==ncol(model.matrix(fit)),formula=paste(deparse(form),collapse=' '))
}
donors <- do.call(rbind,donor_rows)
donors$donor_status <- factor(donors$donor_status,levels=c('ND','Aab+','T1D'))
saveRDS(donor_models,file.path(out,'separate_donor_linear_models.rds'))
write.csv(donors,file.path(out,'donor_coefficients.csv'),row.names=FALSE)
kw <- kruskal.test(standardized_beta~donor_status,data=donors)
write.csv(data.frame(test='Kruskal-Wallis',statistic=unname(kw$statistic),
  df=unname(kw$parameter),p_value=kw$p.value,donors=nrow(donors)),
  file.path(out,'donor_disease_global_test.csv'),row.names=FALSE)
comparisons <- list()
for (pair in combn(c('ND','Aab+','T1D'),2,simplify=FALSE)) {
  x <- donors$standardized_beta[donors$donor_status==pair[1]]
  y <- donors$standardized_beta[donors$donor_status==pair[2]]
  exact_available <- anyDuplicated(c(x,y))==0L
  test <- wilcox.test(x,y,paired=FALSE,alternative='two.sided',exact=exact_available,correct=TRUE)
  comparisons[[length(comparisons)+1]] <- data.frame(group_1=pair[1],group_2=pair[2],
    donors_1=length(x),donors_2=length(y),statistic=unname(test$statistic),p_value=test$p.value,
    method=if(exact_available)'Exact two-sided Wilcoxon rank-sum' else 'Tie-corrected normal approximation with continuity correction',
    median_difference_group2_minus_group1=median(y)-median(x))
}
comparison <- do.call(rbind,comparisons)
comparison$p_bh <- p.adjust(comparison$p_value,method='BH')
write.csv(comparison,file.path(out,'donor_pairwise_disease_tests.csv'),row.names=FALSE)
group_summary <- do.call(rbind,lapply(c('ND','Aab+','T1D'),function(g) {
  z <- donors[donors$donor_status==g,]
  data.frame(group=g,donors=nrow(z),objects=sum(z$objects),mean_beta=mean(z$standardized_beta),
    median_beta=median(z$standardized_beta),q25=quantile(z$standardized_beta,.25),
    q75=quantile(z$standardized_beta,.75))
}))
write.csv(group_summary,file.path(out,'donor_group_summary.csv'),row.names=FALSE)
print(kw); print(comparison); print(group_summary)
write.csv(do.call(rbind,diagnostics),file.path(out,'model_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,groups),file.path(out,'group_slopes.csv'),row.names=FALSE)
write.csv(do.call(rbind,contrasts),file.path(out,'pairwise_slope_differences.csv'),row.names=FALSE)
write.csv(do.call(rbind,interactions),file.path(out,'disease_interaction_tests.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
