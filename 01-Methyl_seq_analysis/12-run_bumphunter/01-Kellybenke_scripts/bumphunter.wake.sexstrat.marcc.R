
library(bsseq)
library(bumphunter)

setwd("/data/NHLBI_BCB/Fayaz/19-rales12-bumphunter-results-dichot")

pheno1 = read.table("rales1.pheno.txt", header=T, row.names=1)
pheno2 = read.table("rales2.pheno.txt", header=T, row.names=1)

out.values=c(pheno1$meanNCortW,pheno2$meanNCortW)
age=c(pheno1$age,pheno2$age)
sex=c(pheno1$sex,pheno2$sex)

boy=which(sex==0)
girl=which(sex==1)

out.boy=out.values[boy]
out.girl=out.values[girl]

age.boy=age[boy]
age.girl=age[girl]

beta = read.table('rales12.betavalues.txt.gz')
print(dim(beta))
beta2=beta[,1:ncol(beta)] 

beta.boy=beta2[,boy]
beta.girl=beta2[,girl]

rownames = rownames(beta)
rownamessplit = unlist(strsplit(rownames, split="_"))
select = grep("chr", rownamessplit)
chr = rownamessplit[select]
pos = rownamessplit[-select]
pos = as.numeric(pos)

#Night
setwd("/data/NHLBI_BCB/Fayaz/19-rales12-bumphunter-results-dichot")
load("rales12.night.sva.Rdata")
sv=as.data.frame(night.sva$sv)
colnames(sv)=paste("sv",1:ncol(sv),sep="")

sv.boy=sv[boy,]
sv.girl=sv[girl,]

#Boy Model

pheno_mat = model.matrix(~as.numeric(out.boy) + as.numeric(age.boy) + as.matrix(sv.boy))
boy.bump = bumphunter(as.matrix(beta.boy),pheno_mat,chr=chr,pos=pos,coef=2,pickCutoff=TRUE,pickCutoffQ=0.99,B=1000,maxGap=300,smooth=FALSE,smoothFunction=loessByCluster,nullMethod="bootstrap")

setwd("/data/NHLBI_BCB/Fayaz/19-rales12-bumphunter-results-dichot")
write.table(boy.bump$table, file="rales12.bump.wake.boy.txt")
boy.bump.pvaluesMarginal = data.frame(cbind(beta=boy.bump$coef, pvaluesMarginal=boy.bump$pvaluesMarginal))
write.table(boy.bump.pvaluesMarginal, file="rales12.bump.pvalues.wake.boy.txt")

#Girl Model

pheno_mat = model.matrix(~as.numeric(out.girl) + as.numeric(age.girl) + as.matrix(sv.girl))
girl.bump = bumphunter(as.matrix(beta.girl),pheno_mat,chr=chr,pos=pos,coef=2,pickCutoff=TRUE,pickCutoffQ=0.99,B=1000,maxGap=300,smooth=FALSE,smoothFunction=loessByCluster,nullMethod="bootstrap")

setwd("/data/NHLBI_BCB/Fayaz/19-rales12-bumphunter-results-dichot")
write.table(girl.bump$table, file="rales12.bump.wake.girl.txt")
girl.bump.pvaluesMarginal = data.frame(cbind(beta=girl.bump$coef, pvaluesMarginal=girl.bump$pvaluesMarginal))
write.table(girl.bump.pvaluesMarginal, file="rales12.bump.pvalues.wake.girl.txt")





