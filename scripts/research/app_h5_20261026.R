out <- 'output/app_hypotheses_20261026/h5'
stopifnot(file.exists(file.path(out,'analysis_specification_before_fit.json')))
d <- read.csv(file.path(out,'data/donor_measurements.csv'))
d$donor_status <- factor(d$donor_status,levels=c('ND','Aab+','T1D'))
om <- list(); pw <- list()
for(key in unique(d$endpoint)){
 z <- d[d$endpoint==key & is.finite(d$value),];fam<-z$family[1]
 allsame <- length(unique(z$value))<2
 test<-if(allsame)NULL else kruskal.test(value~donor_status,data=z)
 om[[key]]<-data.frame(endpoint=key,family=fam,statistic=if(allsame)0 else unname(test$statistic),df=2,p_value=if(allsame)1 else test$p.value,n_donors=nrow(z),status=if(allsame)'all measurements tied' else 'computed')
 for(pp in combn(levels(d$donor_status),2,simplify=FALSE)){
  a<-z$value[z$donor_status==pp[1]];b<-z$value[z$donor_status==pp[2]];tied<-anyDuplicated(c(a,b))>0
  allsame<-length(unique(c(a,b)))<2
  w<-if(allsame)NULL else wilcox.test(a,b,exact=!tied,correct=TRUE)
  pw[[paste(key,pp[1],pp[2])]]<-data.frame(endpoint=key,family=fam,group_1=pp[1],group_2=pp[2],n_1=length(a),n_2=length(b),statistic=if(allsame)NA_real_ else unname(w$statistic),p_value=if(allsame)1 else w$p.value,test=if(allsame)'all measurements tied' else w$method)
 }
}
om<-do.call(rbind,om);pw<-do.call(rbind,pw)
om$p_bh_family<-NA_real_;pw$p_bh_family<-NA_real_
for(fam in c('primary','secondary')){
 i<-om$family==fam;om$p_bh_family[i]<-p.adjust(om$p_value[i],'BH')
 i<-pw$family==fam;pw$p_bh_family[i]<-p.adjust(pw$p_value[i],'BH')
}
write.csv(om,file.path(out,'omnibus_tests.csv'),row.names=FALSE)
write.csv(pw,file.path(out,'pairwise_tests.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
cat('Saved',nrow(om),'omnibus tests and',nrow(pw),'pairwise tests.\n')
