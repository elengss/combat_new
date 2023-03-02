library(lme4)

load("/well/combat/users/nif917/ref/CD8/covid_flu_temra.RData") 

list()->res1
list()->res2
list()->ano

pvec<-vector()
evec<-vector()
for(i in 1:ncol(matnum)){
data.frame(as.numeric(as.character(cd8small[zzz,])),matnum[,i],obssmall["age",],obssmall["pool",],obssmall["id",],as.factor(obssmall["sex",]))->df
colnames(df)<-c("gene","vgene","age","pool","id","sex")
lmer(gene~vgene+age+sex+(1|pool)+(1|id),data=df)->res1[[i]]
lmer(gene~age+sex+(1|pool)+(1|id),data=df)->res2[[i]]
anova(res2[[i]],res1[[i]])->ano[[i]]}

pv<-vector()
for(i in 1:length(ano)){
ano[[i]][["Pr(>Chisq)"]][[2]]->pv[i]
print(i)}

save(pv,file=paste("/well/combat/users/nif917/ref/CD8/results_flu_beta_temra/zzz",".Rd",sep=""))


