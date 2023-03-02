
read.csv("cd8_covid_tem.csv")->cd
read.csv("/well/combat/users/nif917/ref/CD8/CD8_covid_tem/obs.csv")->obs

obs[which(obs$TCR_chain_composition=="double_alpha_beta"),]->tmp
tmp[which(tmp$TCR_doublet==0),]->tmp2
tmp2[which(tmp2$TCR_productive_TRB==1),]->obs2

read.table("/well/combat/users/nif917/ref/out")->out
cbind(obs2$X,obs2$TCR_v_gene_TRB)->x
x[which(is.na(x[,2])==F),]->x2
sapply(strsplit(x2[,2],"\\*"),"[[",1)->x2[,2]
unique(x2[,2])->ug
matrix(nrow=length(ug),ncol=nrow(x2))->mat
rownames(mat)<-ug
for(j in 1:nrow(x2)){
mat[which(ug==x2[j,2]),j]<-"0/1:1.000"}
mat->mat2
mat2[is.na(mat)==T]<-"0/0:0.000"
colnames(mat2)<-x2[,1]
out[grep("TRBV",out$V9),]->out2

c("chr7","142553725","142554208","ENSG00000282040",";","ENSG00000282040",";","+","TRBV7-8",";")->vec1
c("chr7","142362570","142363134","ENSG00000282543",";","ENSG00000282543",";","+","TRBV4-3",";")->vec2
c("chr7","142561449","142562408","ENSG00000282054",";","ENSG00000282054",";","+","TRBV5-8",";")->vec3
c("chr7","142549110","142549542","ENSG00000282610",";","ENSG00000282610",";","+","TRBV6-9",";")->vec4
rbind(out2,vec1,vec2,vec3,vec4)->out3
cbind(out3[,1],out3[,2],out3[,9],"A","G","100","PASS",".","GT:DS")->front

intersect(colnames(mat2),cd[,1])->intercells
mat2[,match(intercells,colnames(mat2))]->mat3
front[match(rownames(mat3),front[,3]),]->front2
cbind(front2,mat3)->frontout
colnames(frontout)[1:9]<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

sapply(strsplit(out[,4],"\\."),"[[",1)->refensg
intersect(refensg,colnames(cd))->inter
cd[,2:ncol(cd)]->cd8new
cd[,1]->rownames(cd8new)
cd8new[,match(inter,colnames(cd8new))]->cd8m
out[match(inter,refensg),]->outm
refensg[match(inter,refensg)]->refm
t(cd8m)->tcd8m
tcd8m[,match(intercells,colnames(tcd8m))]->cd8small

cbind(outm[,1],outm[,2],1+as.numeric(outm[,2]),refm,refm,outm[,8],cd8small)->bed
colnames(bed)[1:6]<-c("Chr","start","end","pid","gid","strand")
###write.table(bed,file="../CD8/pheno_cd8_small.txt",quote=F,col.names=T,row.names=F,sep="\t")

obs[match(intercells,obs[,1]),]->obsm
rep("M",nrow(obsm))->sex
sex[which(obsm$Sex==1)]<-"F"
### convert Sex
rbind(obsm$Age,obsm$Source,obsm$Pool_ID,obsm$COMBAT_ID,sex)->obssmall
rownames(obssmall)<-c("age","dis","pool","id","sex")
colnames(obssmall)<-obsm$X

mat3->matn
matn[which(mat3=="0/0:0.000")]<-0
matn[which(mat3=="0/1:1.000")]<-1
matn->matnum
apply(matn,1,as.numeric)->matnum
prcomp(matnum,scale=T,center=T)->res

#### convert mat to numerical matrix
mat3->matn
matn[which(mat3=="0/0:0.000")]<-0
matn[which(mat3=="0/1:1.000")]<-1
matn->matnum
apply(matn,1,as.numeric)->matnum
prcomp(matnum,scale=T,center=T)->res

t(cd8small)->tcd
prcomp(tcd,scale=T,center=T)->resg

save.image("covid_tem.RData")

