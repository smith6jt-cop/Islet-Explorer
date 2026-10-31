#!/usr/bin/env Rscript
# Original H8: disease-by-continuous-field comparison, donors biological units.
options(warn=1)
out <- 'output/h5_h10_original_plan_20261026/h8'
t <- read.csv(file.path(out,'all_tissue_tiles.csv'),stringsAsFactors=FALSE)
specs <- data.frame(model=c('primary_200um','grid100','grid400','identity_lower','identity_upper','wholecell_size','area_size'),
 grid_um=c(200,100,400,200,200,200,200),
 endpoint=c(rep('small_beta_cluster_observed',3),'small_beta_cluster_lower','small_beta_cluster_upper','small_beta_cluster_wholecells','small_beta_cluster_area'))
write.csv(specs,file.path(out,'model_specifications.csv'),row.names=FALSE)
within_z <- function(v) if(sd(v)>1e-12) (v-mean(v))/sd(v) else rep(0,length(v))
fields <- c('immune','epithelial','stromal')
coefs <- list(); tests <- list(); pairs <- list();fits <- list();scaling <- list()
for (i in seq_len(nrow(specs))) {
 s<-specs[i,];d<-t[t$grid_um==s$grid_um,]
 d$response <- d[[paste0(s$endpoint,'_density_mm2')]]
 d$response_z <- as.numeric(scale(d$response))
 d$immune <- log1p(d$immune_cells_density_mm2)
 d$epithelial <- log1p(d$epithelial_signal_density_mm2)
 d$stromal <- log1p(d$stromal_signal_density_mm2)
 d$area <- log(d$valid_area_um2)
 d$missing_fraction <- d$field_marker_missing_cells/pmax(d$extraendocrine_cells,1)
 scaling[[s$model]] <- data.frame(model=s$model,global_response_mean=mean(d$response),global_response_sd=sd(d$response))
 rec<-list()
 for (case in unique(d$case_id)) {
   z<-d[d$case_id==case,]
   for (v in c(fields,'area','missing_fraction')) z[[paste0(v,'_z')]]<-within_z(z[[v]])
   vars<-paste0(c(fields,'area','missing_fraction'),'_z')
   vars<-vars[sapply(z[vars],sd)>1e-12]
   m<-lm(reformulate(vars,response='response_z'),z)
   if(m$rank<ncol(model.matrix(m)))stop(paste('Rank deficient donor field model',s$model,case))
   fits[[paste(s$model,case,sep='_')]]<-m
   for (field in fields) {
     b<-coef(m)[[paste0(field,'_z')]]
     if(all(z$response==0))b<-0
     rec[[paste(case,field,sep='_')]]<-data.frame(model=s$model,grid_um=s$grid_um,endpoint=s$endpoint,
       case_id=case,donor_status=z$donor_status[1],field=field,standardized_beta=b,
       tiles=nrow(z),zero_object_tiles=sum(z$all_geometric_objects==0),
       zero_endpoint_tiles=sum(z[[s$endpoint]]==0),endpoint_count=sum(z[[s$endpoint]]),
       valid_area_mm2=sum(z$valid_area_um2)/1e6,design_rank=m$rank,design_columns=ncol(model.matrix(m)))
   }
 }
 c<-do.call(rbind,rec);coefs[[s$model]]<-c
 for(field in fields) {
   z<-c[c$field==field,];z$donor_status<-factor(z$donor_status,levels=c('ND','Aab+','T1D'))
   k<-kruskal.test(standardized_beta~donor_status,z)
   tests[[paste(s$model,field)]]<-data.frame(model=s$model,endpoint=s$endpoint,grid_um=s$grid_um,field=field,
     statistic=unname(k$statistic),df=unname(k$parameter),p=k$p.value,donors=nrow(z))
   pr<-list()
   for(gs in combn(levels(z$donor_status),2,simplify=FALSE)) {
     a<-z$standardized_beta[z$donor_status==gs[1]];b<-z$standardized_beta[z$donor_status==gs[2]]
     w<-wilcox.test(a,b,exact=!anyDuplicated(c(a,b)),alternative='two.sided')
     pr[[paste(gs,collapse='_')]]<-data.frame(model=s$model,endpoint=s$endpoint,grid_um=s$grid_um,field=field,
       group1=gs[1],group2=gs[2],n1=length(a),n2=length(b),median1=median(a),median2=median(b),
       p=w$p.value,statistic=unname(w$statistic),method=w$method)
   }
   pr<-do.call(rbind,pr);pr$p_BH<-p.adjust(pr$p,'BH');pairs[[paste(s$model,field)]]<-pr
 }
}
ct<-do.call(rbind,coefs);ot<-do.call(rbind,tests);pt<-do.call(rbind,pairs)
ot$p_BH_three_fields<-ave(ot$p,ot$model,FUN=function(v)p.adjust(v,'BH'))
pt$p_BH_all_nine<-ave(pt$p,pt$model,FUN=function(v)p.adjust(v,'BH'))
write.csv(ct,file.path(out,'donor_field_coefficients.csv'),row.names=FALSE)
write.csv(ot,file.path(out,'disease_field_overall_tests.csv'),row.names=FALSE)
write.csv(pt,file.path(out,'disease_field_pairwise_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,scaling),file.path(out,'response_standardization.csv'),row.names=FALSE)
saveRDS(fits,file.path(out,'donor_field_models.rds'))

# Natural-density comparisons by object class; each point is a whole donor.
census<-read.csv(file.path(out,'donor_census.csv'))
density_rows<-list();density_tests<-list();density_pairs<-list()
for(endpoint in c('small_beta_cluster_observed','small_beta_islet','large_beta_islet','beta_single')) {
 z<-data.frame(case_id=census$case_id,donor_status=factor(census$donor_status,levels=c('ND','Aab+','T1D')),
    endpoint=endpoint,count=census[[endpoint]],area_mm2=census$valid_tissue_area_mm2)
 z$density_mm2<-z$count/z$area_mm2;density_rows[[endpoint]]<-z
 k<-kruskal.test(density_mm2~donor_status,z)
 density_tests[[endpoint]]<-data.frame(endpoint=endpoint,statistic=unname(k$statistic),df=unname(k$parameter),p=k$p.value)
 pr<-list()
 for(gs in combn(levels(z$donor_status),2,simplify=FALSE)) {
   a<-z$density_mm2[z$donor_status==gs[1]];b<-z$density_mm2[z$donor_status==gs[2]]
   w<-wilcox.test(a,b,exact=!anyDuplicated(c(a,b)),alternative='two.sided')
   pr[[paste(gs,collapse='_')]]<-data.frame(endpoint=endpoint,group1=gs[1],group2=gs[2],
     median1=median(a),median2=median(b),p=w$p.value,method=w$method)
 }
 pr<-do.call(rbind,pr);pr$p_BH<-p.adjust(pr$p,'BH');density_pairs[[endpoint]]<-pr
}
write.csv(do.call(rbind,density_rows),file.path(out,'donor_object_class_densities.csv'),row.names=FALSE)
write.csv(do.call(rbind,density_tests),file.path(out,'object_class_density_overall_tests.csv'),row.names=FALSE)
dc<-do.call(rbind,density_pairs);dc$p_BH_all_twelve<-p.adjust(dc$p,'BH')
write.csv(dc,file.path(out,'object_class_density_pairwise_tests.csv'),row.names=FALSE)
print(ot[ot$model=='primary_200um',]);print(pt[pt$model=='primary_200um',])
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
