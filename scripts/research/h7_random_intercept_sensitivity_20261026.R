#!/usr/bin/env Rscript
# Supporting standard mixed-model sensitivity after zero random-slope variance.
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(jsonlite))
out <- 'output/h5_h10_original_plan_20261026/h7'
d <- read.csv(file.path(out,'primary_20um_model_objects.csv'))
d$case_id <- factor(d$case_id)
d$donor_status <- factor(d$donor_status,levels=c('ND','Aab+','T1D'))
form <- response_z ~ exposure_z * donor_status + nerve_area_fraction_z + cd31_area_fraction_z + log_core_area_z + log_perimeter_z + log_nearby_immune_density_z + (1|case_id)
write_json(list(role='Supporting sensitivity only; does not replace declared donor-fallback primary.',
 trigger='The primary independent donor random-slope variance was zero.',
 assumption='No remaining variation of the adjusted proximity slope between donors within a disease state; the reported standard errors use this assumption.',
 formula=paste(deparse(form),collapse=' '),
 endpoint='Exactly the same 946 objects, transformations, exposure and five fixed covariates as primary_20um',
 inference='lmerTest/Satterthwaite, with BH across 3 conditional disease slopes and separately across 3 disease-slope contrasts',
 local_precedent='Barlow2025 and Bodenmiller2025 donor random-intercept mixed regressions'),
 file.path(out,'random_intercept_sensitivity_plan.json'),pretty=TRUE,auto_unbox=TRUE)
fit <- lmer(form,d,REML=TRUE,control=lmerControl(optimizer='bobyqa'))
saveRDS(fit,file.path(out,'primary_20um_random_intercept_sensitivity.rds'))
tr <- emtrends(fit,~donor_status,var='exposure_z',lmer.df='satterthwaite')
ss <- as.data.frame(summary(tr,infer=c(TRUE,TRUE),adjust='none'));ss$p_BH<-p.adjust(ss$p.value,'BH');ss$model<-'random_intercept_20um'
pp <- as.data.frame(summary(pairs(tr),infer=c(TRUE,TRUE),adjust='none'));pp$p_BH<-p.adjust(pp$p.value,'BH');pp$model<-'random_intercept_20um'
aa <- as.data.frame(anova(fit,ddf='Satterthwaite'));aa$term<-rownames(aa)
fx <- as.data.frame(coef(summary(fit)));fx$term<-rownames(fx)
write.csv(ss,file.path(out,'random_intercept_conditional_slopes.csv'),row.names=FALSE)
write.csv(pp,file.path(out,'random_intercept_disease_comparisons.csv'),row.names=FALSE)
write.csv(aa,file.path(out,'random_intercept_model_tests.csv'),row.names=FALSE)
write.csv(fx,file.path(out,'random_intercept_fixed_coefficients.csv'),row.names=FALSE)
write_json(list(singular=isSingular(fit),convergence=fit@optinfo$conv$lme4$messages,variance=as.data.frame(VarCorr(fit))),
 file.path(out,'random_intercept_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
print(ss);print(pp)
