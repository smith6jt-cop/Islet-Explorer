#!/usr/bin/env Rscript
# H1 with statistical methods represented in the user's local paper collection.
suppressPackageStartupMessages({library(lme4);library(lmerTest);library(emmeans)})
out <- 'output/h1_h4_literature_methods_20261025/h1'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(file.path(out,'analysis_objects.csv'),stringsAsFactors=FALSE)
d$case_id <- factor(d$case_id)
d$donor_status <- factor(d$donor_status,levels=c('ND','Aab+','T1D'))
d$object_class <- factor(d$object_class)
stopifnot(nrow(d)==2937,!('78' %in% d$case_id),length(unique(d$case_id))==20)
for(n in c('log_area','log_area_ratio','log_nearby','log_field')) d[[n]] <- as.numeric(scale(d[[n]]))
d$beta_present <- as.integer(d$n_beta_core>0)
d$beta_density_k <- d$beta_density_mm2/1000
control <- glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=150000))
rows <- list(); diagnostics <- list(); models <- list()

fit_one <- function(kind,scope,exclusive=FALSE) {
  z <- d
  if(scope=='T1D_separate') z <- droplevels(subset(z,donor_status=='T1D'))
  if(kind=='positive_density') z <- droplevels(subset(z,n_beta_core>0))
  if(exclusive) {
    z$n_immune_nearby <- z$n_immune_nearby_exclusive
    z <- subset(z,n_immune_nearby>0)
    z$log_nearby <- as.numeric(scale(log1p(z$n_immune_nearby)))
  }
  z$n_immune_peri_model <- z$n_immune_nearby-z$n_immune_core
  term <- if(kind=='reciprocal') 'beta_density_log2' else 'allocation10'
  exposure <- if(scope=='group_interaction') paste0(term,'*donor_status') else
    if(scope=='pooled') paste(term,'+ donor_status') else term
  rhs <- paste(exposure,'+ log_area + log_area_ratio + log_nearby + log_field + object_class + (1|case_id)')
  response <- switch(kind,density='beta_density_mm2',count='n_beta_core',
    presence='beta_present',positive_density='beta_density_k',
    reciprocal='cbind(n_immune_core,n_immune_peri_model)')
  # The freely estimated log_area coefficient already adjusts count for area.
  # Absorbing a log(area) offset into that coefficient is an equivalent design
  # parameterization and avoids the NB optimizer's large-intercept failure.
  offset <- switch(kind,count='',
    reciprocal='+ offset(log(core_area_um2/peri_area_um2))','')
  formula <- as.formula(paste(response,'~',rhs,offset))
  name <- paste(kind,scope,if(exclusive)'exclusive' else 'shared',sep='_')
  caught <- character()
  fit <- tryCatch(withCallingHandlers({
    if(kind=='density') lmerTest::lmer(formula,data=z,REML=TRUE,
      control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=150000)))
    else if(kind=='count') lme4::glmer.nb(formula,data=z,control=control)
    else if(kind=='positive_density') glmer(formula,data=z,family=Gamma(link='log'),nAGQ=9,control=control)
    else glmer(formula,data=z,family=binomial(),control=control)
  },warning=function(w){caught <<- c(caught,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error')) {
    diagnostics[[name]] <<- data.frame(model=name,status='failed',objects=nrow(z),
      donors=length(unique(z$case_id)),singular=NA,warnings=conditionMessage(fit),formula=paste(deparse(formula),collapse=' '))
    print(paste(name,'FAILED',conditionMessage(fit)));return(invisible(NULL))
  }
  saveRDS(fit,file.path(out,paste0(name,'.rds')))
  trend <- emtrends(fit,if(scope=='group_interaction')~donor_status else ~1,var=term,
                   lmer.df='satterthwaite')
  tab <- as.data.frame(summary(trend,infer=c(TRUE,TRUE),adjust='none'))
  estimate <- tab[[paste0(term,'.trend')]]
  lo <- tab[[if('lower.CL'%in%names(tab))'lower.CL' else 'asymp.LCL']]
  hi <- tab[[if('upper.CL'%in%names(tab))'upper.CL' else 'asymp.UCL']]
  status <- if(scope=='group_interaction')as.character(tab$donor_status) else if(scope=='T1D_separate')'T1D' else 'Pooled'
  row <- data.frame(model=name,kind=kind,scope=scope,group=status,estimate=estimate,
    SE=tab$SE,df=tab$df,lower=lo,upper=hi,p=tab$p.value,
    effect=if(kind=='density')estimate else exp(estimate),
    effect_lower=if(kind=='density')lo else exp(lo),
    effect_upper=if(kind=='density')hi else exp(hi),
    objects=sapply(status,function(g)if(g=='Pooled')nrow(z) else sum(z$donor_status==g)),
    donors=sapply(status,function(g)length(unique(z$case_id[if(g=='Pooled')rep(TRUE,nrow(z)) else z$donor_status==g]))))
  rows[[name]] <<- row
  diagnostics[[name]] <<- data.frame(model=name,status='fit',objects=nrow(z),
    donors=length(unique(z$case_id)),singular=isSingular(fit),
    warnings=paste(unique(c(caught,fit@optinfo$conv$lme4$messages)),collapse=' | '),
    formula=paste(deparse(formula),collapse=' '))
  if(scope=='group_interaction') {
    cmp <- as.data.frame(summary(pairs(trend,adjust='BH'),infer=c(TRUE,TRUE)))
    write.csv(cmp,file.path(out,paste0(name,'_group_comparisons.csv')),row.names=FALSE)
    if(kind=='density')write.csv(as.data.frame(anova(fit)),file.path(out,'density_anova.csv'))
  }
  print(row[,c('kind','scope','group','effect','effect_lower','effect_upper','p')]);flush.console()
}

# Declared before estimates: density mixed model is primary; count and
# component models answer specified related questions, never select by p value.
for(kind in c('density','count','presence','positive_density','reciprocal')) {
  for(scope in c('pooled','group_interaction','T1D_separate')) fit_one(kind,scope)
}
fit_one('reciprocal','T1D_separate',exclusive=TRUE)
write.csv(do.call(rbind,rows),file.path(out,'model_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'model_diagnostics.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
