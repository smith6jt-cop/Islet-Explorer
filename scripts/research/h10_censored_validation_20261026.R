#!/usr/bin/env Rscript
# Independent verification of the approved censored-normal likelihood using survival.
suppressPackageStartupMessages(library(survival))
out <- 'output/h5_h10_original_plan_20261026/h10'
d <- read.csv(file.path(out,'donor_inputs.csv'))
d$C_peptide_censored <- tolower(as.character(d$C_peptide_censored))=='true'
d$donor_status <- factor(d$donor_status,levels=c('ND','Aab+','T1D'))
for (v in c('tissue_beta_density_mm2','age_years_npod','admission_glucose_mg_dL','mean_boundary_coverage','mean_immune_localization'))
  d[[paste0(v,'_z')]] <- as.numeric(scale(d[[v]]))
d$spatial_score <- .5*(d$mean_boundary_coverage_z-d$mean_immune_localization_z)
y <- Surv(log(d$C_peptide_ng_mL),!d$C_peptide_censored,type='left')
rows <- list(); coefrows <- list()
for (conditioned in c(FALSE,TRUE)) for (augmented in c(FALSE,TRUE)) {
  label <- if(conditioned)'censored_disease_conditioned' else 'censored'
  name <- if(augmented)'augmented' else 'baseline'
  rhs <- c('tissue_beta_density_mm2_z','age_years_npod_z','admission_glucose_mg_dL_z',
           if(conditioned)'donor_status',if(augmented)'spatial_score')
  fit <- survreg(reformulate(rhs,response='y'),data=d,dist='gaussian',
                control=survreg.control(maxiter=1000,rel.tolerance=1e-12))
  rows[[paste(label,name)]] <- data.frame(analysis=label,model=name,log_likelihood=fit$loglik[2],
    residual_sigma=fit$scale,iterations=fit$iter)
  coefrows[[paste(label,name)]] <- data.frame(analysis=label,model=name,term=names(coef(fit)),estimate=coef(fit))
}
observed <- do.call(rbind,rows)
python <- read.csv(file.path(out,'censored_model_diagnostics.csv'))
check <- merge(observed,python[,c('analysis','model','log_likelihood','residual_sigma')],
               by=c('analysis','model'),suffixes=c('_R','_Python'))
check$log_likelihood_difference <- check$log_likelihood_R-check$log_likelihood_Python
check$sigma_difference <- check$residual_sigma_R-check$residual_sigma_Python
stopifnot(max(abs(check$log_likelihood_difference))<1e-7,max(abs(check$sigma_difference))<1e-6)
write.csv(check,file.path(out,'censored_independent_validation.csv'),row.names=FALSE)
write.csv(do.call(rbind,coefrows),file.path(out,'censored_R_validation_coefficients.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'censored_R_validation_session.txt'))
print(check)
