
# from http://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

setwd("/data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-v4-all-combined/02-extreme_tertiles/meanNCortB/")
phenodata = read.table(file="pheno.txt", header=T, row.names=1)

res = read.csv(file="stat_results_gene_counts_0_1_meanNCortB.csv", header=T, row.names=1)

library(gage)
library(gageData)
library("AnnotationDbi")

# mouse or human
organism_name<-"human"


if(organism_name == "mouse" ){
  library("org.Mm.eg.db")
  columns(org.Mm.eg.db)
  data(go.sets.mm)
  data(go.subs.mm)
}
if(organism_name == "human" ){
  library("org.Hs.eg.db")
  columns(org.Hs.eg.db)
  data(go.sets.hs)
  data(go.subs.hs)
}

library(dplyr)
library(DESeq2)
library(topGO)
library(genefilter)
library(fdrtool)
library(xlsx)

# Gene_COUNT - uses genesymbol for mapping
if(organism_name == "human" ){
  res$symbol = mapIds(org.Hs.eg.db, keys=as.character(res$genesymbol), column="SYMBOL",   keytype="ALIAS", multiVals="first")
  res$ensmbl = mapIds(org.Hs.eg.db, keys=as.character(res$genesymbol), column="ENSEMBL",  keytype="ALIAS", multiVals="first")
  res$entrez = mapIds(org.Hs.eg.db, keys=as.character(res$genesymbol), column="ENTREZID", keytype="ALIAS", multiVals="first")
  res$name =   mapIds(org.Hs.eg.db, keys=as.character(res$genesymbol), column="GENENAME", keytype="ALIAS", multiVals="first")
}
if(organism_name == "mouse" ){
  res$symbol = mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="SYMBOL",   keytype="ALIAS", multiVals="first")
  res$ensmbl = mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="ENSEMBL",  keytype="ALIAS", multiVals="first")
  res$entrez = mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="ENTREZID", keytype="ALIAS", multiVals="first")
  res$name =   mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="GENENAME", keytype="ALIAS", multiVals="first")
}

foldchanges = res$log2FC
names(foldchanges) = res$entrez
head(foldchanges)

if(organism_name == "mouse" ){
  gobpsets = go.sets.mm[go.subs.mm$BP]
}
if(organism_name == "human" ){
  gobpsets = go.sets.hs[go.subs.hs$BP]
}


gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
head(gobpres$greater[, 2:5])
head(gobpres$less[, 2:5])
go_header=c("GO term","stat.mean","p.val","q.val","set.size")
write.table(t(go_header), file="GO-BP.up.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(t(go_header), file="GO-BP.down.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(na.omit(gobpres$greater[, 2:5]), file = "GO-BP.up.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )
write.table(na.omit(gobpres$less[, 2:5]),    file = "GO-BP.down.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )

if(organism_name == "mouse" ){
  gobpsets = go.sets.mm[go.subs.mm$MF]
}
if(organism_name == "human" ){
  gobpsets = go.sets.hs[go.subs.hs$MF]
}


gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
head(gobpres$greater[, 2:5])
head(gobpres$less[, 2:5])
go_header=c("GO term","stat.mean","p.val","q.val","set.size")
write.table(t(go_header), file="GO-MF.up.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(t(go_header), file="GO-MF.down.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(na.omit(gobpres$greater[, 2:5]), file = "GO-MF.up.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )
write.table(na.omit(gobpres$less[, 2:5]),    file = "GO-MF.down.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )

if(organism_name == "mouse" ){
  gobpsets = go.sets.mm[go.subs.mm$CC]
}
if(organism_name == "human" ){
  gobpsets = go.sets.hs[go.subs.hs$CC]
}

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
head(gobpres$greater[, 2:5])
head(gobpres$less[, 2:5])
go_header=c("GO term","stat.mean","p.val","q.val","set.size")
write.table(t(go_header), file="GO-CC.up.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(t(go_header), file="GO-CC.down.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(na.omit(gobpres$greater[, 2:5]), file = "GO-CC.up.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )
write.table(na.omit(gobpres$less[, 2:5]),    file = "GO-CC.down.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )

# get the gene names for top 10
#go=go.gsets(species = "mouse", id.type = "EG")
#go.gs=go$go.sets
#head(go.gs)
#countData <- read.table("10.5_14.5_p_mm10_coding_cleaned_entrez.tsv", header=T, sep="\t", row.names=1)
#for (gs in rownames(gobpres$greater)[1:10]) {
#  outname = gsub(" |:|/", "_", substr(gs, 12, 100))
#  geneData(genes = go.gs, exprs = countData, ref = NULL, 
#           samp = NULL, outname = outname, txt = T, heatmap = T,
#           limit = 3, scatterplot = T)
#  }

# http://genomespot.blogspot.com.au/2014/09/data-analysis-step-8-pathway-analysis.html

xlsx::write.xlsx(read.table("GO-BP.up.sig.tsv", header=T, sep="\t"),file='GO_UP.xlsx',row.names=F,sheetName = 'BP.UP')
xlsx::write.xlsx(read.table("GO-MF.up.sig.tsv", header=T, sep="\t"),file='GO_UP.xlsx',row.names=F,sheetName = 'MF.UP', append=T)
xlsx::write.xlsx(read.table("GO-CC.up.sig.tsv", header=T, sep="\t"),file='GO_UP.xlsx',row.names=F,sheetName = 'CC.UP', append=T)

xlsx::write.xlsx(read.table("GO-BP.down.sig.tsv", header=T, sep="\t"),file='GO_DOWN.xlsx',row.names=F,sheetName = 'BP.DOWN', append=T)
xlsx::write.xlsx(read.table("GO-MF.down.sig.tsv", header=T, sep="\t"),file='GO_DOWN.xlsx',row.names=F,sheetName = 'MF.DOWN', append=T)
xlsx::write.xlsx(read.table("GO-CC.down.sig.tsv", header=T, sep="\t"),file='GO_DOWN.xlsx',row.names=F,sheetName = 'CC.DOWN', append=T)

if(fasle) {
## Heat Map for top results based on criteria sent by Kim
	res_subset_hmap = subset(res, (log2fc < -0.50 | log2fc > 0.50) & pval < 0.05)
  	res_subset_hmap = res_subset_hmap[order(abs(as.numeric(res_subset_hmap$log2fc))),]

	#use top 20% of genes (upregulated and downregulated) if the number of genes in pathway is greater than or equal to 50
	if ((dim(res_subset_hmap)[1]) >= 50) {
		res_subset_hmap_updown_top25 = res_subset_hmap[(dim(res_subset_hmap)[1]-(0.20 * dim(res_subset_hmap)[1])):dim(res_subset_hmap)[1],]
	} else {
		res_subset_hmap_updown_top25 = res_subset_hmap	
	}
	res_subset_hmap_updown_top25 = res_subset_hmap_updown_top25[order(as.numeric(res_subset_hmap_updown_top25$fc)),]

	scaleforheatmap = function(x){
		indexes = which(x==0)
		x[indexes] = mean(x)
		newx = log2(x/mean(x))
		meannewx = mean(newx)
		newxx = newx - meannewx
	}

	xmat = subset(res_subset_hmap_updown_top25, select = grep("S", names(res_subset_hmap)))
  	rownames(xmat) = res_subset_hmap_updown_top25$geneName
  	x = apply(xmat, 1, function(x) scaleforheatmap(as.numeric(x)))	
	newxmat = cbind(t(log2(cpm(x)+1)))
	colnames(newxmat) = c("S10_KO","S6_KO","S7_KO","S8_KO","S9_KO","S1_WT","S2_WT","S3_WT","S4_WT","S5_WT")
	pheatmap(newxmat,filename = paste("pheatmap","_","top_genes_FPKM_pvalue_0.05_log2fc_2updown_transcript_level",".pdf",sep=""),fontsize_row = 6,fontsize_col = 6,fontsize = 5,show_rownames = T, scale="row", cluster_rows=F, cluster_cols=F, col=colorRampPalette( c("green", "black", "red"), space="rgb")(64), cellwidth=10, cellheight=10, main='Top_Genes_FPKM_transcript_level')
#	pheatmap(newxmat,filename = paste("pheatmap","_",pathway,".pdf",sep=""),fontsize_row = 6,fontsize_col = 6,scale="row", cluster_rows=F, cluster_cols=F, col=colorRampPalette( c("green", "black", "red"), space="rgb")(64), cellwidth=10, cellheight=10, fontsize = 4.3,show_rownames = T, main=pathway)
	dev.off()

}



