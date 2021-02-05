
library(bsseq)
library(bumphunter)

setwd("/data/NHLBI_BCB/Fayaz/18-rales12-bumphunter-results")

pheno1 = read.table("rales1.pheno.txt", header=T, row.names=1)
pheno2 = read.table("rales2.pheno.txt", header=T, row.names=1)

out=c(pheno1$meanNCortB,pheno2$meanNCortB)
age=c(pheno1$age,pheno2$age)
sex=c(pheno1$sex,pheno2$sex)


beta = read.table('rales12.betavalues.txt.gz')
print(dim(beta))
beta2=beta[,1:ncol(beta)] 

rownames = rownames(beta)
rownamessplit = unlist(strsplit(rownames, split="_"))
select = grep("chr", rownamessplit)
chr = rownamessplit[select]
pos = rownamessplit[-select]
pos = as.numeric(pos)
#chr=beta[,"chromosome"]
#pos=beta[,"startpos"]

#Night
setwd("/data/NHLBI_BCB/Fayaz/18-rales12-bumphunter-results")
load("rales12.night.sva.Rdata")
sv=as.data.frame(night.sva$sv)
colnames(sv)=paste("sv",1:ncol(sv),sep="")


pheno_mat = model.matrix(~as.numeric(out) + as.numeric(age) + as.factor(sex) + as.matrix(sv) )
r12.bump = bumphunter(as.matrix(beta2),pheno_mat,chr=chr,pos=pos,coef=2,pickCutoff=TRUE,pickCutoffQ=0.99,B=1000,maxGap=300,smooth=FALSE,smoothFunction=loessByCluster,nullMethod="bootstrap")

setwd("/data/NHLBI_BCB/Fayaz/18-rales12-bumphunter-results")
write.table(r12.bump$table, file="rales12.bump.night.txt")
rales12.bump.pvaluesMarginal = data.frame(cbind(beta=r12.bump$coef, pvaluesMarginal=r12.bump$pvaluesMarginal))
write.table(rales12.bump.pvaluesMarginal, file="rales12.bump.pvalues.night.txt")


