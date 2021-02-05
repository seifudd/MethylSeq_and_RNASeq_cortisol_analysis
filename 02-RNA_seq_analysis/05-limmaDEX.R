
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(vsn)
library(devtools)
library(rafalib)
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(plyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library("sva")
library(limma)
library(edgeR)
library(qvalue)

setwd("/data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-v4-all-combined/02-extreme_tertiles/meanNCortB/")
counts_data = read.table(file="counts.txt", header=T, row.names=1, check.names=F)
phenodata = read.table(file="pheno.txt", header=T, row.names=1)

group1 = 0 
group2 = 1
pheno="meanNCortB"
feature = "gene"

phenodata = subset(phenodata, phenodata$group %in% c(group1,group2))
phenodata = phenodata[order(phenodata$group),]
counts_data = subset(counts_data, select=rownames(phenodata))

## filtering
#Keep genes with least 1 count-per-million reads (cpm) in at least 17 samples
isexpr = rowSums(cpm(counts_data)>1) >= 17
table(isexpr)
counts_data = counts_data[isexpr,]
head(cpm(counts_data))

##### sva
## Set null and alternative models (ignore batch)
mod1 = model.matrix(~group, data=phenodata)
mod0 = cbind(mod1[,1])
svseq = svaseq(cpm(counts_data),mod1,mod0)$sv
dim(svseq)

##### voom
dge = DGEList(counts = counts_data)
dge = calcNormFactors(dge)
design = model.matrix(~group + round + sex + Age + svseq[,1] + svseq[,2] + svseq[,3], data=phenodata)
v = voom(dge, design, plot = FALSE, save.plot = TRUE)
fit = lmFit(v, design)
fit = eBayes(fit)
log2FC = fit$coefficients[, 2]
pvalue = fit$p.value[, 2]
coefs = fit$coefficients[, 2]
qvalue = qvalue(pvalue)$qvalues
gene_res = cbind(data.frame(counts_data,log2FC, pvalue, qvalue, check.names=F))


## Significant at FDR 5 and 10%
#results$sig <- results$qvalue < 0.05
#results$sig10 <- results$qvalue < 0.1

## Re-order by qvalue
results = gene_res[order(gene_res$qvalue, decreasing = FALSE), ]
head(results)
annot = read.table(file="../../gencode.v25.primary_assembly.annotation.txt", header=1, row.names=1)
resultsannot = merge(x=annot, y=results, by="row.names")
resultsannot = resultsannot[order(resultsannot$qvalue, decreasing = FALSE), ]
write.csv(resultsannot, paste("stat_results_",feature,"_counts_",group1,"_",group2,"_",pheno,".csv", sep=""), row.names=FALSE)

# pca plot - using results with p-value < 0.05 or any other subset based on significant DE for e.g. fc, log2fc
statres = subset(results, subset=results$pvalue < 0.05)
statres = subset(statres, select = grep("12", names(statres)))

pca = prcomp(t(log2(cpm(statres) +1)), scale=F)
pdf(paste("pca_plot_p_less_than_0.05_pairs_plot_first_five_PCs_",feature,"_",group1,"_",group2,"_",pheno,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",group1,"_",group2,"_",pheno, sep=""), pch=21, bg=c('blue','orange')[as.vector(unclass(as.factor(as.vector(phenodata$group))))], cex.main=0.75)
dev.off()

write.csv(pca$x[,1:5], paste("pca",".csv", sep=""), row.names=TRUE)


# pdf(paste("pca_plot_p_less_than_0.05_pc1_vs_pc2_",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
# plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('blue','orange')[as.vector(unclass(as.factor(as.vector(phenodata$group))))], cex=1.20, cex.main=0.80, xlab="PC1", ylab="PC2")
# text(pca$x[,1], pca$x[,2], labels=phenodata$SampleName, pos=3, cex=0.65)
# legend("topright", inset=c(-0.2,0), legend=c(group1,group2), pch=c(16,16), col=levels(as.factor(c('blue','orange')[as.vector(unclass(as.factor(as.vector(phenodata$group))))])), title="Group", cex=0.70)
# dev.off()

## Plotting abundance distribution, before and after applying filter, transcripts, genes, fpkm
counts_comparison = subset(results, select = grep("12", names(results)))
counts = counts_comparison
colnames(counts) = rownames(phenodata)
#counts = cpm(counts)
#counts = log2(cpm(counts) +1)
counts = log2(counts+1)
pdf(file=paste("log2_cpm_distribution_",feature,"_",group1,"_",group2,"_",pheno,".pdf", sep=""))
boxplot(counts, col=c('blue','orange')[as.vector(unclass(as.factor(as.vector(phenodata$tertile))))], las=2, ylab='log2_cpm', main=paste("log2_cpm distribution:",feature,"_",group1,"_",group2,"_",pheno), xlab="sample", cex.main=0.80, cex.axis=0.50)
dev.off()

## Plotting volcano plot, Differential Expression Analysis-qval
data = results
q.trans = -1*log(data$qvalue, base = 10)
logfold = data$log2FC
pdf(file=paste("volcano_qval_",feature,"_",group1,"_",group2,"_",pheno,".pdf", sep=""))
plot(logfold,q.trans,type="n",ylab="-1*log10(q-value)",xlab="log2(fold change)",main=paste("COUNT:",feature,"_",group1,"_",group2,"_",pheno), cex.main=0.80, xlim=range(logfold), ylim=range(q.trans))
points(logfold,q.trans,col="black",cex=0.65)
points(logfold[(q.trans>1.3&logfold>1.0)],q.trans[(q.trans>1.3&logfold>1.0)],col="red",pch=16,cex=0.65)
points(logfold[(q.trans>1.3&logfold<(-1.0))],q.trans[(q.trans>1.3&logfold<(-1.0))],col="green",pch=16,cex=0.65)
# text(logfold[(q.trans>1.3&logfold>1.0)],q.trans[(q.trans>1.3&logfold>1.0)],labels=as.character(data$geneName[(q.trans>1.3&logfold>1.0)]), cex=0.65)
# text(logfold[(q.trans>1.3&logfold<(-1.0))],q.trans[(q.trans>1.3&logfold<(-1.0))],labels=as.character(data$geneName[(q.trans>1.3&logfold<(-1.0))]), cex=0.65)
abline(h=1.3)
abline(v=-1.0)
abline(v=1.0)
dev.off()

## Plotting volcano plot, Differential Expression Analysis-pval, transcripts, genes, fpkm
p.trans = -1*log(data$pvalue, base = 10)
logfold = data$log2FC
pdf(file=paste("volcano_pval_",feature,"_",group1,"_",group2,"_",pheno,".pdf", sep=""))
plot(logfold,p.trans,type="n",ylab="-1*log10(p-value)",xlab="log2(fold change)",main=paste("COUNT:",feature,"_",group1,"_",group2,"_",pheno), cex.main=0.80, xlim=range(logfold), ylim=range(p.trans))
points(logfold,p.trans,col="black",cex=0.65)
points(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],col="red",pch=16,cex=0.65)
points(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],col="green",pch=16,cex=0.65)
# text(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],labels=as.character(data$geneName[(p.trans>1.3&logfold>1.0)]), cex=0.65)
# text(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],labels=as.character(data$geneName[(p.trans>1.3&logfold<(-1.0))]), cex=0.65)
abline(h=1.3)
abline(v=-1.0)
abline(v=1.0)
dev.off()

pdf(file=paste("qval_distribution_",feature,"_",group1,"_",group2,"_",pheno,".pdf", sep=""))
hist(data$qvalue, col="grey", border="white", xlab="qval", ylab="", main=paste("COUNT:",feature,"_",group1,"_",group2,"_",pheno))
dev.off()
pdf(file=paste("pval_distribution_",feature,"_",group1,"_",group2,"_",pheno,".pdf", sep=""))
hist(data$pvalue, col="grey", border="white", xlab="pval", ylab="", main=paste("COUNT:",feature,"_",group1,"_",group2,"_",pheno))
dev.off()


