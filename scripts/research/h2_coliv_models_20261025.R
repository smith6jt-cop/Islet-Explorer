#!/usr/bin/env Rscript
options(warn=1)
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(emmeans))
out <- 'output/h1_h4_literature_methods_20261025/h2'
d <- read.csv(file.path(out, 'h2_analysis_objects.csv'))
d$case_id <- factor(d$case_id)
d$donor_status <- factor(d$donor_status, levels=c('ND','Aab+','T1D'))
d$object_class <- factor(d$object_class)
d$core_density <- d$n_immune_core * 1e6 / d$core_area_um2
d$log_density <- log1p(d$core_density)
d$peri_density_log <- log1p(d$n_immune_peri * 1e6 / d$peri_area_um2)
d$log_core_area <- log(d$core_area_um2)
d$log_perimeter <- log(d$core_perimeter_um)
d$log_edge_distance <- log1p(pmax(d$core_distance_to_tissue_edge_um, 0))
for (col in c('peri_density_log','log_core_area','log_perimeter','log_edge_distance','capillary_core_area_fraction')) {
  d[[col]] <- as.numeric(scale(d[[col]]))
}
models <- list(); results <- list(); diagnostics <- list()
controls <- 'peri_density_log + log_core_area + log_perimeter + log_edge_distance + capillary_core_area_fraction + object_class'
for (width in c(2,5)) {
  d$coverage_10pp <- d[[paste0('coliv_coverage_',width,'um')]] / .1
  for (outcome in c('core_density','log_density')) {
    for (group in c('Pooled','ND','Aab+','T1D')) {
      z <- if (group=='Pooled') d else d[d$donor_status==group,]
      z <- droplevels(z[is.finite(z$coverage_10pp),])
      form <- as.formula(paste(outcome,'~ coverage_10pp +',controls,
                              if (group=='Pooled') '+ donor_status' else '', '+ (1|case_id)'))
      name <- paste(width,outcome,group,sep='_')
      m <- lmer(form,z,REML=TRUE,control=lmerControl(optimizer='bobyqa'))
      a <- coef(summary(m))['coverage_10pp',]
      q <- qt(.975,a[['df']])
      results[[name]] <- data.frame(width_um=width,outcome=outcome,model='separate_group',group=group,
                        n=nobs(m),donors=length(unique(z$case_id)),estimate=a[['Estimate']],se=a[['Std. Error']],
                        df=a[['df']],lower=a[['Estimate']]-q*a[['Std. Error']],upper=a[['Estimate']]+q*a[['Std. Error']],
                        p=a[['Pr(>|t|)']],singular=isSingular(m), stringsAsFactors=FALSE)
      diagnostics[[name]] <- list(formula=deparse(form),sigma=sigma(m),variance=as.data.frame(VarCorr(m)),
                                  convergence=m@optinfo$conv, excluded_missing_coverage=nrow(if(group=='Pooled')d else d[d$donor_status==group,])-nobs(m))
      models[[name]] <- m
      print(results[[name]])
    }
    form <- as.formula(paste(outcome,'~ coverage_10pp * donor_status +',controls,'+ (1|case_id)'))
    z <- d[is.finite(d$coverage_10pp),]
    m <- lmer(form,z,REML=TRUE,control=lmerControl(optimizer='bobyqa'))
    trend <- as.data.frame(summary(emtrends(m,~donor_status,var='coverage_10pp',lmer.df='satterthwaite'),infer=c(TRUE,TRUE)))
    name <- paste(width,outcome,'interaction',sep='_')
    results[[name]] <- data.frame(width_um=width,outcome=outcome,model='interaction_common_adjustment',group=trend$donor_status,
                                  n=nobs(m),donors=length(unique(z$case_id)),estimate=trend$coverage_10pp.trend,se=trend$SE,
                                  df=trend$df,lower=trend$lower.CL,upper=trend$upper.CL,p=trend$p.value,singular=isSingular(m))
    models[[name]] <- m
    write.csv(as.data.frame(anova(m)),file.path(out,paste0('anova_',name,'.csv')))
  }
}
write.csv(do.call(rbind,results),file.path(out,'h2_model_results.csv'),row.names=FALSE)
saveRDS(models,file.path(out,'h2_models.rds'))
saveRDS(diagnostics,file.path(out,'h2_model_diagnostics.rds'))
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
