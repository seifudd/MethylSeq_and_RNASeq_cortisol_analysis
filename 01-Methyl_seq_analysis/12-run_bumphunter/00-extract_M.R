library(bsseq)
library(bumphunter)

pheno = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/emotional_stress_high_low.txt", header=T, row.names=1)

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
write.table(RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr1/RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=T)
rm(RL_chr1_filtered_10.32.cov)
rm(RL_chr1_filtered_10.32.cpgs)
rm(RL_chr1_filtered_10.32.cpgs.M.df)
rm(RL_chr1_filtered_10.32.cpgs.granges)
rm(RL_chr1_filtered_10.32.cpgs.granges.df)
rm(RL_chr1_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr2/RL_chr2_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr2_filtered_10.32.cov)
rm(RL_chr2_filtered_10.32.cpgs)
rm(RL_chr2_filtered_10.32.cpgs.M.df)
rm(RL_chr2_filtered_10.32.cpgs.granges)
rm(RL_chr2_filtered_10.32.cpgs.granges.df)
rm(RL_chr2_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr3/RL_chr3_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr3_filtered_10.32.cov)
rm(RL_chr3_filtered_10.32.cpgs)
rm(RL_chr3_filtered_10.32.cpgs.M.df)
rm(RL_chr3_filtered_10.32.cpgs.granges)
rm(RL_chr3_filtered_10.32.cpgs.granges.df)
rm(RL_chr3_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr4/RL_chr4_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr4_filtered_10.32.cov)
rm(RL_chr4_filtered_10.32.cpgs)
rm(RL_chr4_filtered_10.32.cpgs.M.df)
rm(RL_chr4_filtered_10.32.cpgs.granges)
rm(RL_chr4_filtered_10.32.cpgs.granges.df)
rm(RL_chr4_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr5/RL_chr5_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr5_filtered_10.32.cov)
rm(RL_chr5_filtered_10.32.cpgs)
rm(RL_chr5_filtered_10.32.cpgs.M.df)
rm(RL_chr5_filtered_10.32.cpgs.granges)
rm(RL_chr5_filtered_10.32.cpgs.granges.df)
rm(RL_chr5_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr6/RL_chr6_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr6_filtered_10.32.cov)
rm(RL_chr6_filtered_10.32.cpgs)
rm(RL_chr6_filtered_10.32.cpgs.M.df)
rm(RL_chr6_filtered_10.32.cpgs.granges)
rm(RL_chr6_filtered_10.32.cpgs.granges.df)
rm(RL_chr6_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr7/RL_chr7_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr7_filtered_10.32.cov)
rm(RL_chr7_filtered_10.32.cpgs)
rm(RL_chr7_filtered_10.32.cpgs.M.df)
rm(RL_chr7_filtered_10.32.cpgs.granges)
rm(RL_chr7_filtered_10.32.cpgs.granges.df)
rm(RL_chr7_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr8/RL_chr8_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr8_filtered_10.32.cov)
rm(RL_chr8_filtered_10.32.cpgs)
rm(RL_chr8_filtered_10.32.cpgs.M.df)
rm(RL_chr8_filtered_10.32.cpgs.granges)
rm(RL_chr8_filtered_10.32.cpgs.granges.df)
rm(RL_chr8_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr9/RL_chr9_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr9_filtered_10.32.cov)
rm(RL_chr9_filtered_10.32.cpgs)
rm(RL_chr9_filtered_10.32.cpgs.M.df)
rm(RL_chr9_filtered_10.32.cpgs.granges)
rm(RL_chr9_filtered_10.32.cpgs.granges.df)
rm(RL_chr9_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr10/RL_chr10_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr10_filtered_10.32.cov)
rm(RL_chr10_filtered_10.32.cpgs)
rm(RL_chr10_filtered_10.32.cpgs.M.df)
rm(RL_chr10_filtered_10.32.cpgs.granges)
rm(RL_chr10_filtered_10.32.cpgs.granges.df)
rm(RL_chr10_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr11/RL_chr11_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr11_filtered_10.32.cov)
rm(RL_chr11_filtered_10.32.cpgs)
rm(RL_chr11_filtered_10.32.cpgs.M.df)
rm(RL_chr11_filtered_10.32.cpgs.granges)
rm(RL_chr11_filtered_10.32.cpgs.granges.df)
rm(RL_chr11_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr12/RL_chr12_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr12_filtered_10.32.cov)
rm(RL_chr12_filtered_10.32.cpgs)
rm(RL_chr12_filtered_10.32.cpgs.M.df)
rm(RL_chr12_filtered_10.32.cpgs.granges)
rm(RL_chr12_filtered_10.32.cpgs.granges.df)
rm(RL_chr12_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr13/RL_chr13_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr13_filtered_10.32.cov)
rm(RL_chr13_filtered_10.32.cpgs)
rm(RL_chr13_filtered_10.32.cpgs.M.df)
rm(RL_chr13_filtered_10.32.cpgs.granges)
rm(RL_chr13_filtered_10.32.cpgs.granges.df)
rm(RL_chr13_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr14/RL_chr14_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr14_filtered_10.32.cov)
rm(RL_chr14_filtered_10.32.cpgs)
rm(RL_chr14_filtered_10.32.cpgs.M.df)
rm(RL_chr14_filtered_10.32.cpgs.granges)
rm(RL_chr14_filtered_10.32.cpgs.granges.df)
rm(RL_chr14_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr15/RL_chr15_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr15_filtered_10.32.cov)
rm(RL_chr15_filtered_10.32.cpgs)
rm(RL_chr15_filtered_10.32.cpgs.M.df)
rm(RL_chr15_filtered_10.32.cpgs.granges)
rm(RL_chr15_filtered_10.32.cpgs.granges.df)
rm(RL_chr15_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr16/RL_chr16_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr16_filtered_10.32.cov)
rm(RL_chr16_filtered_10.32.cpgs)
rm(RL_chr16_filtered_10.32.cpgs.M.df)
rm(RL_chr16_filtered_10.32.cpgs.granges)
rm(RL_chr16_filtered_10.32.cpgs.granges.df)
rm(RL_chr16_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr17/RL_chr17_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr17_filtered_10.32.cov)
rm(RL_chr17_filtered_10.32.cpgs)
rm(RL_chr17_filtered_10.32.cpgs.M.df)
rm(RL_chr17_filtered_10.32.cpgs.granges)
rm(RL_chr17_filtered_10.32.cpgs.granges.df)
rm(RL_chr17_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr18/RL_chr18_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr18_filtered_10.32.cov)
rm(RL_chr18_filtered_10.32.cpgs)
rm(RL_chr18_filtered_10.32.cpgs.M.df)
rm(RL_chr18_filtered_10.32.cpgs.granges)
rm(RL_chr18_filtered_10.32.cpgs.granges.df)
rm(RL_chr18_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr19/RL_chr19_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr19_filtered_10.32.cov)
rm(RL_chr19_filtered_10.32.cpgs)
rm(RL_chr19_filtered_10.32.cpgs.M.df)
rm(RL_chr19_filtered_10.32.cpgs.granges)
rm(RL_chr19_filtered_10.32.cpgs.granges.df)
rm(RL_chr19_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr20/RL_chr20_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr20_filtered_10.32.cov)
rm(RL_chr20_filtered_10.32.cpgs)
rm(RL_chr20_filtered_10.32.cpgs.M.df)
rm(RL_chr20_filtered_10.32.cpgs.granges)
rm(RL_chr20_filtered_10.32.cpgs.granges.df)
rm(RL_chr20_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr21/RL_chr21_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr21_filtered_10.32.cov)
rm(RL_chr21_filtered_10.32.cpgs)
rm(RL_chr21_filtered_10.32.cpgs.M.df)
rm(RL_chr21_filtered_10.32.cpgs.granges)
rm(RL_chr21_filtered_10.32.cpgs.granges.df)
rm(RL_chr21_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chr22/RL_chr22_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chr22_filtered_10.32.cov)
rm(RL_chr22_filtered_10.32.cpgs)
rm(RL_chr22_filtered_10.32.cpgs.M.df)
rm(RL_chr22_filtered_10.32.cpgs.granges)
rm(RL_chr22_filtered_10.32.cpgs.granges.df)
rm(RL_chr22_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chrX/RL_chrX_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chrX_filtered_10.32.cov)
rm(RL_chrX_filtered_10.32.cpgs)
rm(RL_chrX_filtered_10.32.cpgs.M.df)
rm(RL_chrX_filtered_10.32.cpgs.granges)
rm(RL_chrX_filtered_10.32.cpgs.granges.df)
rm(RL_chrX_filtered_10.32.cpgs.M.granges.df)

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
write.table(RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df[,5:dim(RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df)[2]], file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/chrY/RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.txt", row.names=F, col.names=F)
rm(RL_chrY_filtered_10.32.cov)
rm(RL_chrY_filtered_10.32.cpgs)
rm(RL_chrY_filtered_10.32.cpgs.M.df)
rm(RL_chrY_filtered_10.32.cpgs.granges)
rm(RL_chrY_filtered_10.32.cpgs.granges.df)
rm(RL_chrY_filtered_10.32.cpgs.M.granges.df)

rm(RL_chrY_filtered_10.32.cpgs.M.granges.without.SNPs.df)

