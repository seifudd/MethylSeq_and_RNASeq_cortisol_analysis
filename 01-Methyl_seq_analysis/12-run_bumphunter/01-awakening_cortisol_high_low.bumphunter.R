library(bsseq)
library(bumphunter)

pheno = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/awakening_cortisol_high_low.txt", header=T, row.names=1)
load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/00-runSVA/RL_awakening_cortisol_filtered_10.32.cpgs.M.granges.without.SNPs.svaobj.Rdata")
pheno_mat = model.matrix(~pheno[,"dx"] + RL_awakening_cortisol_filtered_10.32.cpgs.M.granges.without.SNPs.svaobj$sv)


load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr1/RL_chr1_filtered_10.32.rda")
col.names = sampleNames(RL_chr1_filtered_10.32)
RL_chr1_filtered_10.32.cov = getCoverage(RL_chr1_filtered_10.32, type="Cov")
colnames(RL_chr1_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr1_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr1_filtered_10.32.cpgs = RL_chr1_filtered_10.32[keepLoci.ex,]
RL_chr1_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr1_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr1_filtered_10.32.cpgs.M.df) = col.names
RL_chr1_filtered_10.32.cpgs.M.df = subset(RL_chr1_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr1_filtered_10.32.cpgs.M.df = apply(RL_chr1_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr1_filtered_10.32.cpgs.M.df = t(RL_chr1_filtered_10.32.cpgs.M.df)
RL_chr1_filtered_10.32.cpgs.granges = granges(RL_chr1_filtered_10.32.cpgs)
RL_chr1_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr1_filtered_10.32.cpgs.granges),startpos=start(RL_chr1_filtered_10.32.cpgs.granges),endpos=end(RL_chr1_filtered_10.32.cpgs.granges),strand=strand(RL_chr1_filtered_10.32.cpgs.granges))
rm(RL_chr1_filtered_10.32)
RL_chr1_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr1_filtered_10.32.cpgs.granges.df, RL_chr1_filtered_10.32.cpgs.M.df), row.names=RL_chr1_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr1 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr1/chr1.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr1 = as.numeric(row.names(RL_chr1_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr1[,1]
RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr1_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr1,] 
rm(RL_chr1_filtered_10.32.cov)
rm(RL_chr1_filtered_10.32.cpgs)
rm(RL_chr1_filtered_10.32.cpgs.M.df)
rm(RL_chr1_filtered_10.32.cpgs.granges)
rm(RL_chr1_filtered_10.32.cpgs.granges.df)
rm(RL_chr1_filtered_10.32.cpgs.M.granges.df)

RL_chr1.bumps = bumphunter(as.matrix(RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr1.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr1/awakening_cortisol_high_low_chr1_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr1_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr1.bumps$coef, coefficients_smooth=RL_chr1.bumps$fitted, pvaluesMarginal=RL_chr1.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr1_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr1/awakening_cortisol_high_low_chr1_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr2/RL_chr2_filtered_10.32.rda")
col.names = sampleNames(RL_chr2_filtered_10.32)
RL_chr2_filtered_10.32.cov = getCoverage(RL_chr2_filtered_10.32, type="Cov")
colnames(RL_chr2_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr2_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr2_filtered_10.32.cpgs = RL_chr2_filtered_10.32[keepLoci.ex,]
RL_chr2_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr2_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr2_filtered_10.32.cpgs.M.df) = col.names
RL_chr2_filtered_10.32.cpgs.M.df = subset(RL_chr2_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr2_filtered_10.32.cpgs.M.df = apply(RL_chr2_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr2_filtered_10.32.cpgs.M.df = t(RL_chr2_filtered_10.32.cpgs.M.df)
RL_chr2_filtered_10.32.cpgs.granges = granges(RL_chr2_filtered_10.32.cpgs)
RL_chr2_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr2_filtered_10.32.cpgs.granges),startpos=start(RL_chr2_filtered_10.32.cpgs.granges),endpos=end(RL_chr2_filtered_10.32.cpgs.granges),strand=strand(RL_chr2_filtered_10.32.cpgs.granges))
rm(RL_chr2_filtered_10.32)
RL_chr2_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr2_filtered_10.32.cpgs.granges.df, RL_chr2_filtered_10.32.cpgs.M.df), row.names=RL_chr2_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr2 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr2/chr2.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr2 = as.numeric(row.names(RL_chr2_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr2[,1]
RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr2_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr2,] 
rm(RL_chr2_filtered_10.32.cov)
rm(RL_chr2_filtered_10.32.cpgs)
rm(RL_chr2_filtered_10.32.cpgs.M.df)
rm(RL_chr2_filtered_10.32.cpgs.granges)
rm(RL_chr2_filtered_10.32.cpgs.granges.df)
rm(RL_chr2_filtered_10.32.cpgs.M.granges.df)

RL_chr2.bumps = bumphunter(as.matrix(RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr2.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr2/awakening_cortisol_high_low_chr2_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr2_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr2.bumps$coef, coefficients_smooth=RL_chr2.bumps$fitted, pvaluesMarginal=RL_chr2.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr2_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr2/awakening_cortisol_high_low_chr2_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr3/RL_chr3_filtered_10.32.rda")
col.names = sampleNames(RL_chr3_filtered_10.32)
RL_chr3_filtered_10.32.cov = getCoverage(RL_chr3_filtered_10.32, type="Cov")
colnames(RL_chr3_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr3_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr3_filtered_10.32.cpgs = RL_chr3_filtered_10.32[keepLoci.ex,]
RL_chr3_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr3_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr3_filtered_10.32.cpgs.M.df) = col.names
RL_chr3_filtered_10.32.cpgs.M.df = subset(RL_chr3_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr3_filtered_10.32.cpgs.M.df = apply(RL_chr3_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr3_filtered_10.32.cpgs.M.df = t(RL_chr3_filtered_10.32.cpgs.M.df)
RL_chr3_filtered_10.32.cpgs.granges = granges(RL_chr3_filtered_10.32.cpgs)
RL_chr3_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr3_filtered_10.32.cpgs.granges),startpos=start(RL_chr3_filtered_10.32.cpgs.granges),endpos=end(RL_chr3_filtered_10.32.cpgs.granges),strand=strand(RL_chr3_filtered_10.32.cpgs.granges))
rm(RL_chr3_filtered_10.32)
RL_chr3_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr3_filtered_10.32.cpgs.granges.df, RL_chr3_filtered_10.32.cpgs.M.df), row.names=RL_chr3_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr3 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr3/chr3.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr3 = as.numeric(row.names(RL_chr3_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr3[,1]
RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr3_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr3,] 
rm(RL_chr3_filtered_10.32.cov)
rm(RL_chr3_filtered_10.32.cpgs)
rm(RL_chr3_filtered_10.32.cpgs.M.df)
rm(RL_chr3_filtered_10.32.cpgs.granges)
rm(RL_chr3_filtered_10.32.cpgs.granges.df)
rm(RL_chr3_filtered_10.32.cpgs.M.granges.df)

RL_chr3.bumps = bumphunter(as.matrix(RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr3.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr3/awakening_cortisol_high_low_chr3_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr3_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr3.bumps$coef, coefficients_smooth=RL_chr3.bumps$fitted, pvaluesMarginal=RL_chr3.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr3_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr3/awakening_cortisol_high_low_chr3_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr4/RL_chr4_filtered_10.32.rda")
col.names = sampleNames(RL_chr4_filtered_10.32)
RL_chr4_filtered_10.32.cov = getCoverage(RL_chr4_filtered_10.32, type="Cov")
colnames(RL_chr4_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr4_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr4_filtered_10.32.cpgs = RL_chr4_filtered_10.32[keepLoci.ex,]
RL_chr4_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr4_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr4_filtered_10.32.cpgs.M.df) = col.names
RL_chr4_filtered_10.32.cpgs.M.df = subset(RL_chr4_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr4_filtered_10.32.cpgs.M.df = apply(RL_chr4_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr4_filtered_10.32.cpgs.M.df = t(RL_chr4_filtered_10.32.cpgs.M.df)
RL_chr4_filtered_10.32.cpgs.granges = granges(RL_chr4_filtered_10.32.cpgs)
RL_chr4_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr4_filtered_10.32.cpgs.granges),startpos=start(RL_chr4_filtered_10.32.cpgs.granges),endpos=end(RL_chr4_filtered_10.32.cpgs.granges),strand=strand(RL_chr4_filtered_10.32.cpgs.granges))
rm(RL_chr4_filtered_10.32)
RL_chr4_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr4_filtered_10.32.cpgs.granges.df, RL_chr4_filtered_10.32.cpgs.M.df), row.names=RL_chr4_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr4 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr4/chr4.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr4 = as.numeric(row.names(RL_chr4_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr4[,1]
RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr4_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr4,] 
rm(RL_chr4_filtered_10.32.cov)
rm(RL_chr4_filtered_10.32.cpgs)
rm(RL_chr4_filtered_10.32.cpgs.M.df)
rm(RL_chr4_filtered_10.32.cpgs.granges)
rm(RL_chr4_filtered_10.32.cpgs.granges.df)
rm(RL_chr4_filtered_10.32.cpgs.M.granges.df)

RL_chr4.bumps = bumphunter(as.matrix(RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr4.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr4/awakening_cortisol_high_low_chr4_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr4_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr4.bumps$coef, coefficients_smooth=RL_chr4.bumps$fitted, pvaluesMarginal=RL_chr4.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr4_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr4/awakening_cortisol_high_low_chr4_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr5/RL_chr5_filtered_10.32.rda")
col.names = sampleNames(RL_chr5_filtered_10.32)
RL_chr5_filtered_10.32.cov = getCoverage(RL_chr5_filtered_10.32, type="Cov")
colnames(RL_chr5_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr5_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr5_filtered_10.32.cpgs = RL_chr5_filtered_10.32[keepLoci.ex,]
RL_chr5_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr5_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr5_filtered_10.32.cpgs.M.df) = col.names
RL_chr5_filtered_10.32.cpgs.M.df = subset(RL_chr5_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr5_filtered_10.32.cpgs.M.df = apply(RL_chr5_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr5_filtered_10.32.cpgs.M.df = t(RL_chr5_filtered_10.32.cpgs.M.df)
RL_chr5_filtered_10.32.cpgs.granges = granges(RL_chr5_filtered_10.32.cpgs)
RL_chr5_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr5_filtered_10.32.cpgs.granges),startpos=start(RL_chr5_filtered_10.32.cpgs.granges),endpos=end(RL_chr5_filtered_10.32.cpgs.granges),strand=strand(RL_chr5_filtered_10.32.cpgs.granges))
rm(RL_chr5_filtered_10.32)
RL_chr5_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr5_filtered_10.32.cpgs.granges.df, RL_chr5_filtered_10.32.cpgs.M.df), row.names=RL_chr5_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr5 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr5/chr5.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr5 = as.numeric(row.names(RL_chr5_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr5[,1]
RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr5_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr5,] 
rm(RL_chr5_filtered_10.32.cov)
rm(RL_chr5_filtered_10.32.cpgs)
rm(RL_chr5_filtered_10.32.cpgs.M.df)
rm(RL_chr5_filtered_10.32.cpgs.granges)
rm(RL_chr5_filtered_10.32.cpgs.granges.df)
rm(RL_chr5_filtered_10.32.cpgs.M.granges.df)

RL_chr5.bumps = bumphunter(as.matrix(RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr5.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr5/awakening_cortisol_high_low_chr5_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr5_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr5.bumps$coef, coefficients_smooth=RL_chr5.bumps$fitted, pvaluesMarginal=RL_chr5.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr5_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr5/awakening_cortisol_high_low_chr5_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr6/RL_chr6_filtered_10.32.rda")
col.names = sampleNames(RL_chr6_filtered_10.32)
RL_chr6_filtered_10.32.cov = getCoverage(RL_chr6_filtered_10.32, type="Cov")
colnames(RL_chr6_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr6_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr6_filtered_10.32.cpgs = RL_chr6_filtered_10.32[keepLoci.ex,]
RL_chr6_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr6_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr6_filtered_10.32.cpgs.M.df) = col.names
RL_chr6_filtered_10.32.cpgs.M.df = subset(RL_chr6_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr6_filtered_10.32.cpgs.M.df = apply(RL_chr6_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr6_filtered_10.32.cpgs.M.df = t(RL_chr6_filtered_10.32.cpgs.M.df)
RL_chr6_filtered_10.32.cpgs.granges = granges(RL_chr6_filtered_10.32.cpgs)
RL_chr6_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr6_filtered_10.32.cpgs.granges),startpos=start(RL_chr6_filtered_10.32.cpgs.granges),endpos=end(RL_chr6_filtered_10.32.cpgs.granges),strand=strand(RL_chr6_filtered_10.32.cpgs.granges))
rm(RL_chr6_filtered_10.32)
RL_chr6_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr6_filtered_10.32.cpgs.granges.df, RL_chr6_filtered_10.32.cpgs.M.df), row.names=RL_chr6_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr6 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr6/chr6.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr6 = as.numeric(row.names(RL_chr6_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr6[,1]
RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr6_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr6,] 
rm(RL_chr6_filtered_10.32.cov)
rm(RL_chr6_filtered_10.32.cpgs)
rm(RL_chr6_filtered_10.32.cpgs.M.df)
rm(RL_chr6_filtered_10.32.cpgs.granges)
rm(RL_chr6_filtered_10.32.cpgs.granges.df)
rm(RL_chr6_filtered_10.32.cpgs.M.granges.df)

RL_chr6.bumps = bumphunter(as.matrix(RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr6.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr6/awakening_cortisol_high_low_chr6_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr6_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr6.bumps$coef, coefficients_smooth=RL_chr6.bumps$fitted, pvaluesMarginal=RL_chr6.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr6_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr6/awakening_cortisol_high_low_chr6_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr7/RL_chr7_filtered_10.32.rda")
col.names = sampleNames(RL_chr7_filtered_10.32)
RL_chr7_filtered_10.32.cov = getCoverage(RL_chr7_filtered_10.32, type="Cov")
colnames(RL_chr7_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr7_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr7_filtered_10.32.cpgs = RL_chr7_filtered_10.32[keepLoci.ex,]
RL_chr7_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr7_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr7_filtered_10.32.cpgs.M.df) = col.names
RL_chr7_filtered_10.32.cpgs.M.df = subset(RL_chr7_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr7_filtered_10.32.cpgs.M.df = apply(RL_chr7_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr7_filtered_10.32.cpgs.M.df = t(RL_chr7_filtered_10.32.cpgs.M.df)
RL_chr7_filtered_10.32.cpgs.granges = granges(RL_chr7_filtered_10.32.cpgs)
RL_chr7_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr7_filtered_10.32.cpgs.granges),startpos=start(RL_chr7_filtered_10.32.cpgs.granges),endpos=end(RL_chr7_filtered_10.32.cpgs.granges),strand=strand(RL_chr7_filtered_10.32.cpgs.granges))
rm(RL_chr7_filtered_10.32)
RL_chr7_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr7_filtered_10.32.cpgs.granges.df, RL_chr7_filtered_10.32.cpgs.M.df), row.names=RL_chr7_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr7 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr7/chr7.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr7 = as.numeric(row.names(RL_chr7_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr7[,1]
RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr7_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr7,] 
rm(RL_chr7_filtered_10.32.cov)
rm(RL_chr7_filtered_10.32.cpgs)
rm(RL_chr7_filtered_10.32.cpgs.M.df)
rm(RL_chr7_filtered_10.32.cpgs.granges)
rm(RL_chr7_filtered_10.32.cpgs.granges.df)
rm(RL_chr7_filtered_10.32.cpgs.M.granges.df)

RL_chr7.bumps = bumphunter(as.matrix(RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr7.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr7/awakening_cortisol_high_low_chr7_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr7_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr7.bumps$coef, coefficients_smooth=RL_chr7.bumps$fitted, pvaluesMarginal=RL_chr7.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr7_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr7/awakening_cortisol_high_low_chr7_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr8/RL_chr8_filtered_10.32.rda")
col.names = sampleNames(RL_chr8_filtered_10.32)
RL_chr8_filtered_10.32.cov = getCoverage(RL_chr8_filtered_10.32, type="Cov")
colnames(RL_chr8_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr8_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr8_filtered_10.32.cpgs = RL_chr8_filtered_10.32[keepLoci.ex,]
RL_chr8_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr8_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr8_filtered_10.32.cpgs.M.df) = col.names
RL_chr8_filtered_10.32.cpgs.M.df = subset(RL_chr8_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr8_filtered_10.32.cpgs.M.df = apply(RL_chr8_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr8_filtered_10.32.cpgs.M.df = t(RL_chr8_filtered_10.32.cpgs.M.df)
RL_chr8_filtered_10.32.cpgs.granges = granges(RL_chr8_filtered_10.32.cpgs)
RL_chr8_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr8_filtered_10.32.cpgs.granges),startpos=start(RL_chr8_filtered_10.32.cpgs.granges),endpos=end(RL_chr8_filtered_10.32.cpgs.granges),strand=strand(RL_chr8_filtered_10.32.cpgs.granges))
rm(RL_chr8_filtered_10.32)
RL_chr8_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr8_filtered_10.32.cpgs.granges.df, RL_chr8_filtered_10.32.cpgs.M.df), row.names=RL_chr8_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr8 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr8/chr8.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr8 = as.numeric(row.names(RL_chr8_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr8[,1]
RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr8_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr8,] 
rm(RL_chr8_filtered_10.32.cov)
rm(RL_chr8_filtered_10.32.cpgs)
rm(RL_chr8_filtered_10.32.cpgs.M.df)
rm(RL_chr8_filtered_10.32.cpgs.granges)
rm(RL_chr8_filtered_10.32.cpgs.granges.df)
rm(RL_chr8_filtered_10.32.cpgs.M.granges.df)

RL_chr8.bumps = bumphunter(as.matrix(RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr8.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr8/awakening_cortisol_high_low_chr8_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr8_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr8.bumps$coef, coefficients_smooth=RL_chr8.bumps$fitted, pvaluesMarginal=RL_chr8.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr8_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr8/awakening_cortisol_high_low_chr8_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr9/RL_chr9_filtered_10.32.rda")
col.names = sampleNames(RL_chr9_filtered_10.32)
RL_chr9_filtered_10.32.cov = getCoverage(RL_chr9_filtered_10.32, type="Cov")
colnames(RL_chr9_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr9_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr9_filtered_10.32.cpgs = RL_chr9_filtered_10.32[keepLoci.ex,]
RL_chr9_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr9_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr9_filtered_10.32.cpgs.M.df) = col.names
RL_chr9_filtered_10.32.cpgs.M.df = subset(RL_chr9_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr9_filtered_10.32.cpgs.M.df = apply(RL_chr9_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr9_filtered_10.32.cpgs.M.df = t(RL_chr9_filtered_10.32.cpgs.M.df)
RL_chr9_filtered_10.32.cpgs.granges = granges(RL_chr9_filtered_10.32.cpgs)
RL_chr9_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr9_filtered_10.32.cpgs.granges),startpos=start(RL_chr9_filtered_10.32.cpgs.granges),endpos=end(RL_chr9_filtered_10.32.cpgs.granges),strand=strand(RL_chr9_filtered_10.32.cpgs.granges))
rm(RL_chr9_filtered_10.32)
RL_chr9_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr9_filtered_10.32.cpgs.granges.df, RL_chr9_filtered_10.32.cpgs.M.df), row.names=RL_chr9_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr9 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr9/chr9.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr9 = as.numeric(row.names(RL_chr9_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr9[,1]
RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr9_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr9,] 
rm(RL_chr9_filtered_10.32.cov)
rm(RL_chr9_filtered_10.32.cpgs)
rm(RL_chr9_filtered_10.32.cpgs.M.df)
rm(RL_chr9_filtered_10.32.cpgs.granges)
rm(RL_chr9_filtered_10.32.cpgs.granges.df)
rm(RL_chr9_filtered_10.32.cpgs.M.granges.df)

RL_chr9.bumps = bumphunter(as.matrix(RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr9.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr9/awakening_cortisol_high_low_chr9_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr9_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr9.bumps$coef, coefficients_smooth=RL_chr9.bumps$fitted, pvaluesMarginal=RL_chr9.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr9_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr9/awakening_cortisol_high_low_chr9_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr10/RL_chr10_filtered_10.32.rda")
col.names = sampleNames(RL_chr10_filtered_10.32)
RL_chr10_filtered_10.32.cov = getCoverage(RL_chr10_filtered_10.32, type="Cov")
colnames(RL_chr10_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr10_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr10_filtered_10.32.cpgs = RL_chr10_filtered_10.32[keepLoci.ex,]
RL_chr10_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr10_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr10_filtered_10.32.cpgs.M.df) = col.names
RL_chr10_filtered_10.32.cpgs.M.df = subset(RL_chr10_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr10_filtered_10.32.cpgs.M.df = apply(RL_chr10_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr10_filtered_10.32.cpgs.M.df = t(RL_chr10_filtered_10.32.cpgs.M.df)
RL_chr10_filtered_10.32.cpgs.granges = granges(RL_chr10_filtered_10.32.cpgs)
RL_chr10_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr10_filtered_10.32.cpgs.granges),startpos=start(RL_chr10_filtered_10.32.cpgs.granges),endpos=end(RL_chr10_filtered_10.32.cpgs.granges),strand=strand(RL_chr10_filtered_10.32.cpgs.granges))
rm(RL_chr10_filtered_10.32)
RL_chr10_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr10_filtered_10.32.cpgs.granges.df, RL_chr10_filtered_10.32.cpgs.M.df), row.names=RL_chr10_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr10 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr10/chr10.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr10 = as.numeric(row.names(RL_chr10_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr10[,1]
RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr10_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr10,] 
rm(RL_chr10_filtered_10.32.cov)
rm(RL_chr10_filtered_10.32.cpgs)
rm(RL_chr10_filtered_10.32.cpgs.M.df)
rm(RL_chr10_filtered_10.32.cpgs.granges)
rm(RL_chr10_filtered_10.32.cpgs.granges.df)
rm(RL_chr10_filtered_10.32.cpgs.M.granges.df)

RL_chr10.bumps = bumphunter(as.matrix(RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr10.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr10/awakening_cortisol_high_low_chr10_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr10_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr10.bumps$coef, coefficients_smooth=RL_chr10.bumps$fitted, pvaluesMarginal=RL_chr10.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr10_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr10/awakening_cortisol_high_low_chr10_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr11/RL_chr11_filtered_10.32.rda")
col.names = sampleNames(RL_chr11_filtered_10.32)
RL_chr11_filtered_10.32.cov = getCoverage(RL_chr11_filtered_10.32, type="Cov")
colnames(RL_chr11_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr11_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr11_filtered_10.32.cpgs = RL_chr11_filtered_10.32[keepLoci.ex,]
RL_chr11_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr11_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr11_filtered_10.32.cpgs.M.df) = col.names
RL_chr11_filtered_10.32.cpgs.M.df = subset(RL_chr11_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr11_filtered_10.32.cpgs.M.df = apply(RL_chr11_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr11_filtered_10.32.cpgs.M.df = t(RL_chr11_filtered_10.32.cpgs.M.df)
RL_chr11_filtered_10.32.cpgs.granges = granges(RL_chr11_filtered_10.32.cpgs)
RL_chr11_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr11_filtered_10.32.cpgs.granges),startpos=start(RL_chr11_filtered_10.32.cpgs.granges),endpos=end(RL_chr11_filtered_10.32.cpgs.granges),strand=strand(RL_chr11_filtered_10.32.cpgs.granges))
rm(RL_chr11_filtered_10.32)
RL_chr11_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr11_filtered_10.32.cpgs.granges.df, RL_chr11_filtered_10.32.cpgs.M.df), row.names=RL_chr11_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr11 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr11/chr11.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr11 = as.numeric(row.names(RL_chr11_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr11[,1]
RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr11_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr11,] 
rm(RL_chr11_filtered_10.32.cov)
rm(RL_chr11_filtered_10.32.cpgs)
rm(RL_chr11_filtered_10.32.cpgs.M.df)
rm(RL_chr11_filtered_10.32.cpgs.granges)
rm(RL_chr11_filtered_10.32.cpgs.granges.df)
rm(RL_chr11_filtered_10.32.cpgs.M.granges.df)

RL_chr11.bumps = bumphunter(as.matrix(RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr11.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr11/awakening_cortisol_high_low_chr11_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr11_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr11.bumps$coef, coefficients_smooth=RL_chr11.bumps$fitted, pvaluesMarginal=RL_chr11.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr11_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr11/awakening_cortisol_high_low_chr11_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr12/RL_chr12_filtered_10.32.rda")
col.names = sampleNames(RL_chr12_filtered_10.32)
RL_chr12_filtered_10.32.cov = getCoverage(RL_chr12_filtered_10.32, type="Cov")
colnames(RL_chr12_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr12_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr12_filtered_10.32.cpgs = RL_chr12_filtered_10.32[keepLoci.ex,]
RL_chr12_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr12_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr12_filtered_10.32.cpgs.M.df) = col.names
RL_chr12_filtered_10.32.cpgs.M.df = subset(RL_chr12_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr12_filtered_10.32.cpgs.M.df = apply(RL_chr12_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr12_filtered_10.32.cpgs.M.df = t(RL_chr12_filtered_10.32.cpgs.M.df)
RL_chr12_filtered_10.32.cpgs.granges = granges(RL_chr12_filtered_10.32.cpgs)
RL_chr12_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr12_filtered_10.32.cpgs.granges),startpos=start(RL_chr12_filtered_10.32.cpgs.granges),endpos=end(RL_chr12_filtered_10.32.cpgs.granges),strand=strand(RL_chr12_filtered_10.32.cpgs.granges))
rm(RL_chr12_filtered_10.32)
RL_chr12_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr12_filtered_10.32.cpgs.granges.df, RL_chr12_filtered_10.32.cpgs.M.df), row.names=RL_chr12_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr12 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr12/chr12.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr12 = as.numeric(row.names(RL_chr12_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr12[,1]
RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr12_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr12,] 
rm(RL_chr12_filtered_10.32.cov)
rm(RL_chr12_filtered_10.32.cpgs)
rm(RL_chr12_filtered_10.32.cpgs.M.df)
rm(RL_chr12_filtered_10.32.cpgs.granges)
rm(RL_chr12_filtered_10.32.cpgs.granges.df)
rm(RL_chr12_filtered_10.32.cpgs.M.granges.df)

RL_chr12.bumps = bumphunter(as.matrix(RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr12.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr12/awakening_cortisol_high_low_chr12_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr12_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr12.bumps$coef, coefficients_smooth=RL_chr12.bumps$fitted, pvaluesMarginal=RL_chr12.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr12_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr12/awakening_cortisol_high_low_chr12_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr13/RL_chr13_filtered_10.32.rda")
col.names = sampleNames(RL_chr13_filtered_10.32)
RL_chr13_filtered_10.32.cov = getCoverage(RL_chr13_filtered_10.32, type="Cov")
colnames(RL_chr13_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr13_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr13_filtered_10.32.cpgs = RL_chr13_filtered_10.32[keepLoci.ex,]
RL_chr13_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr13_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr13_filtered_10.32.cpgs.M.df) = col.names
RL_chr13_filtered_10.32.cpgs.M.df = subset(RL_chr13_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr13_filtered_10.32.cpgs.M.df = apply(RL_chr13_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr13_filtered_10.32.cpgs.M.df = t(RL_chr13_filtered_10.32.cpgs.M.df)
RL_chr13_filtered_10.32.cpgs.granges = granges(RL_chr13_filtered_10.32.cpgs)
RL_chr13_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr13_filtered_10.32.cpgs.granges),startpos=start(RL_chr13_filtered_10.32.cpgs.granges),endpos=end(RL_chr13_filtered_10.32.cpgs.granges),strand=strand(RL_chr13_filtered_10.32.cpgs.granges))
rm(RL_chr13_filtered_10.32)
RL_chr13_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr13_filtered_10.32.cpgs.granges.df, RL_chr13_filtered_10.32.cpgs.M.df), row.names=RL_chr13_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr13 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr13/chr13.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr13 = as.numeric(row.names(RL_chr13_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr13[,1]
RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr13_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr13,] 
rm(RL_chr13_filtered_10.32.cov)
rm(RL_chr13_filtered_10.32.cpgs)
rm(RL_chr13_filtered_10.32.cpgs.M.df)
rm(RL_chr13_filtered_10.32.cpgs.granges)
rm(RL_chr13_filtered_10.32.cpgs.granges.df)
rm(RL_chr13_filtered_10.32.cpgs.M.granges.df)

RL_chr13.bumps = bumphunter(as.matrix(RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr13.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr13/awakening_cortisol_high_low_chr13_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr13_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr13.bumps$coef, coefficients_smooth=RL_chr13.bumps$fitted, pvaluesMarginal=RL_chr13.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr13_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr13/awakening_cortisol_high_low_chr13_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr14/RL_chr14_filtered_10.32.rda")
col.names = sampleNames(RL_chr14_filtered_10.32)
RL_chr14_filtered_10.32.cov = getCoverage(RL_chr14_filtered_10.32, type="Cov")
colnames(RL_chr14_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr14_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr14_filtered_10.32.cpgs = RL_chr14_filtered_10.32[keepLoci.ex,]
RL_chr14_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr14_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr14_filtered_10.32.cpgs.M.df) = col.names
RL_chr14_filtered_10.32.cpgs.M.df = subset(RL_chr14_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr14_filtered_10.32.cpgs.M.df = apply(RL_chr14_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr14_filtered_10.32.cpgs.M.df = t(RL_chr14_filtered_10.32.cpgs.M.df)
RL_chr14_filtered_10.32.cpgs.granges = granges(RL_chr14_filtered_10.32.cpgs)
RL_chr14_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr14_filtered_10.32.cpgs.granges),startpos=start(RL_chr14_filtered_10.32.cpgs.granges),endpos=end(RL_chr14_filtered_10.32.cpgs.granges),strand=strand(RL_chr14_filtered_10.32.cpgs.granges))
rm(RL_chr14_filtered_10.32)
RL_chr14_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr14_filtered_10.32.cpgs.granges.df, RL_chr14_filtered_10.32.cpgs.M.df), row.names=RL_chr14_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr14 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr14/chr14.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr14 = as.numeric(row.names(RL_chr14_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr14[,1]
RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr14_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr14,] 
rm(RL_chr14_filtered_10.32.cov)
rm(RL_chr14_filtered_10.32.cpgs)
rm(RL_chr14_filtered_10.32.cpgs.M.df)
rm(RL_chr14_filtered_10.32.cpgs.granges)
rm(RL_chr14_filtered_10.32.cpgs.granges.df)
rm(RL_chr14_filtered_10.32.cpgs.M.granges.df)

RL_chr14.bumps = bumphunter(as.matrix(RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr14.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr14/awakening_cortisol_high_low_chr14_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr14_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr14.bumps$coef, coefficients_smooth=RL_chr14.bumps$fitted, pvaluesMarginal=RL_chr14.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr14_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr14/awakening_cortisol_high_low_chr14_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr15/RL_chr15_filtered_10.32.rda")
col.names = sampleNames(RL_chr15_filtered_10.32)
RL_chr15_filtered_10.32.cov = getCoverage(RL_chr15_filtered_10.32, type="Cov")
colnames(RL_chr15_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr15_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr15_filtered_10.32.cpgs = RL_chr15_filtered_10.32[keepLoci.ex,]
RL_chr15_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr15_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr15_filtered_10.32.cpgs.M.df) = col.names
RL_chr15_filtered_10.32.cpgs.M.df = subset(RL_chr15_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr15_filtered_10.32.cpgs.M.df = apply(RL_chr15_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr15_filtered_10.32.cpgs.M.df = t(RL_chr15_filtered_10.32.cpgs.M.df)
RL_chr15_filtered_10.32.cpgs.granges = granges(RL_chr15_filtered_10.32.cpgs)
RL_chr15_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr15_filtered_10.32.cpgs.granges),startpos=start(RL_chr15_filtered_10.32.cpgs.granges),endpos=end(RL_chr15_filtered_10.32.cpgs.granges),strand=strand(RL_chr15_filtered_10.32.cpgs.granges))
rm(RL_chr15_filtered_10.32)
RL_chr15_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr15_filtered_10.32.cpgs.granges.df, RL_chr15_filtered_10.32.cpgs.M.df), row.names=RL_chr15_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr15 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr15/chr15.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr15 = as.numeric(row.names(RL_chr15_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr15[,1]
RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr15_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr15,] 
rm(RL_chr15_filtered_10.32.cov)
rm(RL_chr15_filtered_10.32.cpgs)
rm(RL_chr15_filtered_10.32.cpgs.M.df)
rm(RL_chr15_filtered_10.32.cpgs.granges)
rm(RL_chr15_filtered_10.32.cpgs.granges.df)
rm(RL_chr15_filtered_10.32.cpgs.M.granges.df)

RL_chr15.bumps = bumphunter(as.matrix(RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr15.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr15/awakening_cortisol_high_low_chr15_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr15_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr15.bumps$coef, coefficients_smooth=RL_chr15.bumps$fitted, pvaluesMarginal=RL_chr15.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr15_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr15/awakening_cortisol_high_low_chr15_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr16/RL_chr16_filtered_10.32.rda")
col.names = sampleNames(RL_chr16_filtered_10.32)
RL_chr16_filtered_10.32.cov = getCoverage(RL_chr16_filtered_10.32, type="Cov")
colnames(RL_chr16_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr16_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr16_filtered_10.32.cpgs = RL_chr16_filtered_10.32[keepLoci.ex,]
RL_chr16_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr16_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr16_filtered_10.32.cpgs.M.df) = col.names
RL_chr16_filtered_10.32.cpgs.M.df = subset(RL_chr16_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr16_filtered_10.32.cpgs.M.df = apply(RL_chr16_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr16_filtered_10.32.cpgs.M.df = t(RL_chr16_filtered_10.32.cpgs.M.df)
RL_chr16_filtered_10.32.cpgs.granges = granges(RL_chr16_filtered_10.32.cpgs)
RL_chr16_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr16_filtered_10.32.cpgs.granges),startpos=start(RL_chr16_filtered_10.32.cpgs.granges),endpos=end(RL_chr16_filtered_10.32.cpgs.granges),strand=strand(RL_chr16_filtered_10.32.cpgs.granges))
rm(RL_chr16_filtered_10.32)
RL_chr16_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr16_filtered_10.32.cpgs.granges.df, RL_chr16_filtered_10.32.cpgs.M.df), row.names=RL_chr16_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr16 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr16/chr16.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr16 = as.numeric(row.names(RL_chr16_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr16[,1]
RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr16_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr16,] 
rm(RL_chr16_filtered_10.32.cov)
rm(RL_chr16_filtered_10.32.cpgs)
rm(RL_chr16_filtered_10.32.cpgs.M.df)
rm(RL_chr16_filtered_10.32.cpgs.granges)
rm(RL_chr16_filtered_10.32.cpgs.granges.df)
rm(RL_chr16_filtered_10.32.cpgs.M.granges.df)

RL_chr16.bumps = bumphunter(as.matrix(RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr16.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr16/awakening_cortisol_high_low_chr16_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr16_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr16.bumps$coef, coefficients_smooth=RL_chr16.bumps$fitted, pvaluesMarginal=RL_chr16.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr16_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr16/awakening_cortisol_high_low_chr16_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr17/RL_chr17_filtered_10.32.rda")
col.names = sampleNames(RL_chr17_filtered_10.32)
RL_chr17_filtered_10.32.cov = getCoverage(RL_chr17_filtered_10.32, type="Cov")
colnames(RL_chr17_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr17_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr17_filtered_10.32.cpgs = RL_chr17_filtered_10.32[keepLoci.ex,]
RL_chr17_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr17_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr17_filtered_10.32.cpgs.M.df) = col.names
RL_chr17_filtered_10.32.cpgs.M.df = subset(RL_chr17_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr17_filtered_10.32.cpgs.M.df = apply(RL_chr17_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr17_filtered_10.32.cpgs.M.df = t(RL_chr17_filtered_10.32.cpgs.M.df)
RL_chr17_filtered_10.32.cpgs.granges = granges(RL_chr17_filtered_10.32.cpgs)
RL_chr17_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr17_filtered_10.32.cpgs.granges),startpos=start(RL_chr17_filtered_10.32.cpgs.granges),endpos=end(RL_chr17_filtered_10.32.cpgs.granges),strand=strand(RL_chr17_filtered_10.32.cpgs.granges))
rm(RL_chr17_filtered_10.32)
RL_chr17_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr17_filtered_10.32.cpgs.granges.df, RL_chr17_filtered_10.32.cpgs.M.df), row.names=RL_chr17_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr17 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr17/chr17.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr17 = as.numeric(row.names(RL_chr17_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr17[,1]
RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr17_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr17,] 
rm(RL_chr17_filtered_10.32.cov)
rm(RL_chr17_filtered_10.32.cpgs)
rm(RL_chr17_filtered_10.32.cpgs.M.df)
rm(RL_chr17_filtered_10.32.cpgs.granges)
rm(RL_chr17_filtered_10.32.cpgs.granges.df)
rm(RL_chr17_filtered_10.32.cpgs.M.granges.df)

RL_chr17.bumps = bumphunter(as.matrix(RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr17.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr17/awakening_cortisol_high_low_chr17_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr17_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr17.bumps$coef, coefficients_smooth=RL_chr17.bumps$fitted, pvaluesMarginal=RL_chr17.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr17_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr17/awakening_cortisol_high_low_chr17_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr18/RL_chr18_filtered_10.32.rda")
col.names = sampleNames(RL_chr18_filtered_10.32)
RL_chr18_filtered_10.32.cov = getCoverage(RL_chr18_filtered_10.32, type="Cov")
colnames(RL_chr18_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr18_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr18_filtered_10.32.cpgs = RL_chr18_filtered_10.32[keepLoci.ex,]
RL_chr18_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr18_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr18_filtered_10.32.cpgs.M.df) = col.names
RL_chr18_filtered_10.32.cpgs.M.df = subset(RL_chr18_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr18_filtered_10.32.cpgs.M.df = apply(RL_chr18_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr18_filtered_10.32.cpgs.M.df = t(RL_chr18_filtered_10.32.cpgs.M.df)
RL_chr18_filtered_10.32.cpgs.granges = granges(RL_chr18_filtered_10.32.cpgs)
RL_chr18_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr18_filtered_10.32.cpgs.granges),startpos=start(RL_chr18_filtered_10.32.cpgs.granges),endpos=end(RL_chr18_filtered_10.32.cpgs.granges),strand=strand(RL_chr18_filtered_10.32.cpgs.granges))
rm(RL_chr18_filtered_10.32)
RL_chr18_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr18_filtered_10.32.cpgs.granges.df, RL_chr18_filtered_10.32.cpgs.M.df), row.names=RL_chr18_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr18 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr18/chr18.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr18 = as.numeric(row.names(RL_chr18_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr18[,1]
RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr18_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr18,] 
rm(RL_chr18_filtered_10.32.cov)
rm(RL_chr18_filtered_10.32.cpgs)
rm(RL_chr18_filtered_10.32.cpgs.M.df)
rm(RL_chr18_filtered_10.32.cpgs.granges)
rm(RL_chr18_filtered_10.32.cpgs.granges.df)
rm(RL_chr18_filtered_10.32.cpgs.M.granges.df)

RL_chr18.bumps = bumphunter(as.matrix(RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr18.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr18/awakening_cortisol_high_low_chr18_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr18_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr18.bumps$coef, coefficients_smooth=RL_chr18.bumps$fitted, pvaluesMarginal=RL_chr18.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr18_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr18/awakening_cortisol_high_low_chr18_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr19/RL_chr19_filtered_10.32.rda")
col.names = sampleNames(RL_chr19_filtered_10.32)
RL_chr19_filtered_10.32.cov = getCoverage(RL_chr19_filtered_10.32, type="Cov")
colnames(RL_chr19_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr19_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr19_filtered_10.32.cpgs = RL_chr19_filtered_10.32[keepLoci.ex,]
RL_chr19_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr19_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr19_filtered_10.32.cpgs.M.df) = col.names
RL_chr19_filtered_10.32.cpgs.M.df = subset(RL_chr19_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr19_filtered_10.32.cpgs.M.df = apply(RL_chr19_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr19_filtered_10.32.cpgs.M.df = t(RL_chr19_filtered_10.32.cpgs.M.df)
RL_chr19_filtered_10.32.cpgs.granges = granges(RL_chr19_filtered_10.32.cpgs)
RL_chr19_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr19_filtered_10.32.cpgs.granges),startpos=start(RL_chr19_filtered_10.32.cpgs.granges),endpos=end(RL_chr19_filtered_10.32.cpgs.granges),strand=strand(RL_chr19_filtered_10.32.cpgs.granges))
rm(RL_chr19_filtered_10.32)
RL_chr19_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr19_filtered_10.32.cpgs.granges.df, RL_chr19_filtered_10.32.cpgs.M.df), row.names=RL_chr19_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr19 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr19/chr19.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr19 = as.numeric(row.names(RL_chr19_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr19[,1]
RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr19_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr19,] 
rm(RL_chr19_filtered_10.32.cov)
rm(RL_chr19_filtered_10.32.cpgs)
rm(RL_chr19_filtered_10.32.cpgs.M.df)
rm(RL_chr19_filtered_10.32.cpgs.granges)
rm(RL_chr19_filtered_10.32.cpgs.granges.df)
rm(RL_chr19_filtered_10.32.cpgs.M.granges.df)

RL_chr19.bumps = bumphunter(as.matrix(RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr19.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr19/awakening_cortisol_high_low_chr19_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr19_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr19.bumps$coef, coefficients_smooth=RL_chr19.bumps$fitted, pvaluesMarginal=RL_chr19.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr19_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr19/awakening_cortisol_high_low_chr19_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr20/RL_chr20_filtered_10.32.rda")
col.names = sampleNames(RL_chr20_filtered_10.32)
RL_chr20_filtered_10.32.cov = getCoverage(RL_chr20_filtered_10.32, type="Cov")
colnames(RL_chr20_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr20_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr20_filtered_10.32.cpgs = RL_chr20_filtered_10.32[keepLoci.ex,]
RL_chr20_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr20_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr20_filtered_10.32.cpgs.M.df) = col.names
RL_chr20_filtered_10.32.cpgs.M.df = subset(RL_chr20_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr20_filtered_10.32.cpgs.M.df = apply(RL_chr20_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr20_filtered_10.32.cpgs.M.df = t(RL_chr20_filtered_10.32.cpgs.M.df)
RL_chr20_filtered_10.32.cpgs.granges = granges(RL_chr20_filtered_10.32.cpgs)
RL_chr20_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr20_filtered_10.32.cpgs.granges),startpos=start(RL_chr20_filtered_10.32.cpgs.granges),endpos=end(RL_chr20_filtered_10.32.cpgs.granges),strand=strand(RL_chr20_filtered_10.32.cpgs.granges))
rm(RL_chr20_filtered_10.32)
RL_chr20_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr20_filtered_10.32.cpgs.granges.df, RL_chr20_filtered_10.32.cpgs.M.df), row.names=RL_chr20_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr20 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr20/chr20.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr20 = as.numeric(row.names(RL_chr20_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr20[,1]
RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr20_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr20,] 
rm(RL_chr20_filtered_10.32.cov)
rm(RL_chr20_filtered_10.32.cpgs)
rm(RL_chr20_filtered_10.32.cpgs.M.df)
rm(RL_chr20_filtered_10.32.cpgs.granges)
rm(RL_chr20_filtered_10.32.cpgs.granges.df)
rm(RL_chr20_filtered_10.32.cpgs.M.granges.df)

RL_chr20.bumps = bumphunter(as.matrix(RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr20.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr20/awakening_cortisol_high_low_chr20_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr20_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr20.bumps$coef, coefficients_smooth=RL_chr20.bumps$fitted, pvaluesMarginal=RL_chr20.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr20_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr20/awakening_cortisol_high_low_chr20_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr21/RL_chr21_filtered_10.32.rda")
col.names = sampleNames(RL_chr21_filtered_10.32)
RL_chr21_filtered_10.32.cov = getCoverage(RL_chr21_filtered_10.32, type="Cov")
colnames(RL_chr21_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr21_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr21_filtered_10.32.cpgs = RL_chr21_filtered_10.32[keepLoci.ex,]
RL_chr21_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr21_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr21_filtered_10.32.cpgs.M.df) = col.names
RL_chr21_filtered_10.32.cpgs.M.df = subset(RL_chr21_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr21_filtered_10.32.cpgs.M.df = apply(RL_chr21_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr21_filtered_10.32.cpgs.M.df = t(RL_chr21_filtered_10.32.cpgs.M.df)
RL_chr21_filtered_10.32.cpgs.granges = granges(RL_chr21_filtered_10.32.cpgs)
RL_chr21_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr21_filtered_10.32.cpgs.granges),startpos=start(RL_chr21_filtered_10.32.cpgs.granges),endpos=end(RL_chr21_filtered_10.32.cpgs.granges),strand=strand(RL_chr21_filtered_10.32.cpgs.granges))
rm(RL_chr21_filtered_10.32)
RL_chr21_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr21_filtered_10.32.cpgs.granges.df, RL_chr21_filtered_10.32.cpgs.M.df), row.names=RL_chr21_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr21 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr21/chr21.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr21 = as.numeric(row.names(RL_chr21_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr21[,1]
RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr21_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr21,] 
rm(RL_chr21_filtered_10.32.cov)
rm(RL_chr21_filtered_10.32.cpgs)
rm(RL_chr21_filtered_10.32.cpgs.M.df)
rm(RL_chr21_filtered_10.32.cpgs.granges)
rm(RL_chr21_filtered_10.32.cpgs.granges.df)
rm(RL_chr21_filtered_10.32.cpgs.M.granges.df)

RL_chr21.bumps = bumphunter(as.matrix(RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr21.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr21/awakening_cortisol_high_low_chr21_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr21_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr21.bumps$coef, coefficients_smooth=RL_chr21.bumps$fitted, pvaluesMarginal=RL_chr21.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr21_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr21/awakening_cortisol_high_low_chr21_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chr22/RL_chr22_filtered_10.32.rda")
col.names = sampleNames(RL_chr22_filtered_10.32)
RL_chr22_filtered_10.32.cov = getCoverage(RL_chr22_filtered_10.32, type="Cov")
colnames(RL_chr22_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chr22_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chr22_filtered_10.32.cpgs = RL_chr22_filtered_10.32[keepLoci.ex,]
RL_chr22_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chr22_filtered_10.32.cpgs, type="raw"))
colnames(RL_chr22_filtered_10.32.cpgs.M.df) = col.names
RL_chr22_filtered_10.32.cpgs.M.df = subset(RL_chr22_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chr22_filtered_10.32.cpgs.M.df = apply(RL_chr22_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chr22_filtered_10.32.cpgs.M.df = t(RL_chr22_filtered_10.32.cpgs.M.df)
RL_chr22_filtered_10.32.cpgs.granges = granges(RL_chr22_filtered_10.32.cpgs)
RL_chr22_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chr22_filtered_10.32.cpgs.granges),startpos=start(RL_chr22_filtered_10.32.cpgs.granges),endpos=end(RL_chr22_filtered_10.32.cpgs.granges),strand=strand(RL_chr22_filtered_10.32.cpgs.granges))
rm(RL_chr22_filtered_10.32)
RL_chr22_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chr22_filtered_10.32.cpgs.granges.df, RL_chr22_filtered_10.32.cpgs.M.df), row.names=RL_chr22_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chr22 = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chr22/chr22.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chr22 = as.numeric(row.names(RL_chr22_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chr22[,1]
RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chr22_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chr22,] 
rm(RL_chr22_filtered_10.32.cov)
rm(RL_chr22_filtered_10.32.cpgs)
rm(RL_chr22_filtered_10.32.cpgs.M.df)
rm(RL_chr22_filtered_10.32.cpgs.granges)
rm(RL_chr22_filtered_10.32.cpgs.granges.df)
rm(RL_chr22_filtered_10.32.cpgs.M.granges.df)

RL_chr22.bumps = bumphunter(as.matrix(RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chr22.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr22/awakening_cortisol_high_low_chr22_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chr22_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chr22.bumps$coef, coefficients_smooth=RL_chr22.bumps$fitted, pvaluesMarginal=RL_chr22.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chr22_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr22/awakening_cortisol_high_low_chr22_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chrX/RL_chrX_filtered_10.32.rda")
col.names = sampleNames(RL_chrX_filtered_10.32)
RL_chrX_filtered_10.32.cov = getCoverage(RL_chrX_filtered_10.32, type="Cov")
colnames(RL_chrX_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chrX_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chrX_filtered_10.32.cpgs = RL_chrX_filtered_10.32[keepLoci.ex,]
RL_chrX_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chrX_filtered_10.32.cpgs, type="raw"))
colnames(RL_chrX_filtered_10.32.cpgs.M.df) = col.names
RL_chrX_filtered_10.32.cpgs.M.df = subset(RL_chrX_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chrX_filtered_10.32.cpgs.M.df = apply(RL_chrX_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chrX_filtered_10.32.cpgs.M.df = t(RL_chrX_filtered_10.32.cpgs.M.df)
RL_chrX_filtered_10.32.cpgs.granges = granges(RL_chrX_filtered_10.32.cpgs)
RL_chrX_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chrX_filtered_10.32.cpgs.granges),startpos=start(RL_chrX_filtered_10.32.cpgs.granges),endpos=end(RL_chrX_filtered_10.32.cpgs.granges),strand=strand(RL_chrX_filtered_10.32.cpgs.granges))
rm(RL_chrX_filtered_10.32)
RL_chrX_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chrX_filtered_10.32.cpgs.granges.df, RL_chrX_filtered_10.32.cpgs.M.df), row.names=RL_chrX_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chrX = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chrX/chrX.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chrX = as.numeric(row.names(RL_chrX_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chrX[,1]
RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chrX_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chrX,] 
rm(RL_chrX_filtered_10.32.cov)
rm(RL_chrX_filtered_10.32.cpgs)
rm(RL_chrX_filtered_10.32.cpgs.M.df)
rm(RL_chrX_filtered_10.32.cpgs.granges)
rm(RL_chrX_filtered_10.32.cpgs.granges.df)
rm(RL_chrX_filtered_10.32.cpgs.M.granges.df)

RL_chrX.bumps = bumphunter(as.matrix(RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chrX.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chrX/awakening_cortisol_high_low_chrX_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chrX_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chrX.bumps$coef, coefficients_smooth=RL_chrX.bumps$fitted, pvaluesMarginal=RL_chrX.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chrX_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chrX/awakening_cortisol_high_low_chrX_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/14-bsseq_R_data_object_all_samples_divided_into_chromosomes_filtered/chrY/RL_chrY_filtered_10.32.rda")
col.names = sampleNames(RL_chrY_filtered_10.32)
RL_chrY_filtered_10.32.cov = getCoverage(RL_chrY_filtered_10.32, type="Cov")
colnames(RL_chrY_filtered_10.32.cov) = col.names
keepLoci.ex = which(rowSums(RL_chrY_filtered_10.32.cov[,] >= 10) >= 32)
length(keepLoci.ex)
RL_chrY_filtered_10.32.cpgs = RL_chrY_filtered_10.32[keepLoci.ex,]
RL_chrY_filtered_10.32.cpgs.M.df = data.frame(getMeth(RL_chrY_filtered_10.32.cpgs, type="raw"))
colnames(RL_chrY_filtered_10.32.cpgs.M.df) = col.names
RL_chrY_filtered_10.32.cpgs.M.df = subset(RL_chrY_filtered_10.32.cpgs.M.df, select=rownames(pheno))
RL_chrY_filtered_10.32.cpgs.M.df = apply(RL_chrY_filtered_10.32.cpgs.M.df, 1, function(x){x[x==0] = (0+0.0000001); x})
RL_chrY_filtered_10.32.cpgs.M.df = t(RL_chrY_filtered_10.32.cpgs.M.df)
RL_chrY_filtered_10.32.cpgs.granges = granges(RL_chrY_filtered_10.32.cpgs)
RL_chrY_filtered_10.32.cpgs.granges.df = data.frame(chromosome=seqnames(RL_chrY_filtered_10.32.cpgs.granges),startpos=start(RL_chrY_filtered_10.32.cpgs.granges),endpos=end(RL_chrY_filtered_10.32.cpgs.granges),strand=strand(RL_chrY_filtered_10.32.cpgs.granges))
rm(RL_chrY_filtered_10.32)
RL_chrY_filtered_10.32.cpgs.M.granges.df = data.frame(cbind(RL_chrY_filtered_10.32.cpgs.granges.df, RL_chrY_filtered_10.32.cpgs.M.df), row.names=RL_chrY_filtered_10.32.cpgs.granges.df$startpos)
cpgs_with_snps_chrY = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/09-gwillour-methylseq-suicide-project/14-dbSNP142_common_SNPs/chrY/chrY.cpg.dbsnp.C.position.G.position.or.both.actual.txt", header=F)
exclude_cpgs_with_snps_chrY = as.numeric(row.names(RL_chrY_filtered_10.32.cpgs.M.granges.df)) %in% cpgs_with_snps_chrY[,1]
RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df = RL_chrY_filtered_10.32.cpgs.M.granges.df[!exclude_cpgs_with_snps_chrY,] 
rm(RL_chrY_filtered_10.32.cov)
rm(RL_chrY_filtered_10.32.cpgs)
rm(RL_chrY_filtered_10.32.cpgs.M.df)
rm(RL_chrY_filtered_10.32.cpgs.granges)
rm(RL_chrY_filtered_10.32.cpgs.granges.df)
rm(RL_chrY_filtered_10.32.cpgs.M.granges.df)

RL_chrY.bumps = bumphunter(as.matrix(RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]]), design=pheno_mat, chr=RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df$chromosome, pos=RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df$startpos, coef=2, pickCutoff=TRUE, pickCutoffQ=0.99, B=1000, maxGap=300, smooth=FALSE, smoothFunction=loessByCluster, nullMethod="bootstrap")
gc()
write.table(RL_chrY.bumps$table, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chrY/awakening_cortisol_high_low_chrY_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
awakening_cortisol_high_low_chrY_bumps_coef_pvaluesMarginal_smoothcoef = data.frame(cbind(coefficients=RL_chrY.bumps$coef, coefficients_smooth=RL_chrY.bumps$fitted, pvaluesMarginal=RL_chrY.bumps$pvaluesMarginal))
write.table(awakening_cortisol_high_low_chrY_bumps_coef_pvaluesMarginal_smoothcoef, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chrY/awakening_cortisol_high_low_chrY_bumps_coef_pvaluesMarginal_smoothcoef_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt")
rm(RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df)

