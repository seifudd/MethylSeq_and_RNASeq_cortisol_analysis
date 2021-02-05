library(bsseq)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL1.rda")
RL1
RL1_cov = getCoverage(RL1, type="Cov")
sum(RL1_cov >= 10)
round(colMeans(RL1_cov), 1)
RL1_keep_loci = which(RL1_cov >= 10)
RL1_10.X = RL1[RL1_keep_loci,]
round(colMeans(getCoverage(RL1_10.X, type="Cov")), 1)
range(getMeth(RL1_10.X, type="raw"))
rm(RL1)
rm(RL1_10.X)
rm(RL1_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL1_on_target.rda")
RL1_on_target
RL1_on_target_cov = getCoverage(RL1_on_target, type="Cov")
sum(RL1_on_target_cov >= 10)
round(colMeans(RL1_on_target_cov), 1)
RL1_on_target_keep_loci = which(RL1_on_target_cov >= 10)
RL1_on_target_10.X = RL1_on_target[RL1_on_target_keep_loci,]
round(colMeans(getCoverage(RL1_on_target_10.X, type="Cov")), 1)
range(getMeth(RL1_on_target_10.X, type="raw"))
rm(RL1_on_target)
rm(RL1_on_target_10.X)
rm(RL1_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL2.rda")
RL2
RL2_cov = getCoverage(RL2, type="Cov")
sum(RL2_cov >= 10)
round(colMeans(RL2_cov), 1)
RL2_keep_loci = which(RL2_cov >= 10)
RL2_10.X = RL2[RL2_keep_loci,]
round(colMeans(getCoverage(RL2_10.X, type="Cov")), 1)
range(getMeth(RL2_10.X, type="raw"))
rm(RL2)
rm(RL2_10.X)
rm(RL2_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL2_on_target.rda")
RL2_on_target
RL2_on_target_cov = getCoverage(RL2_on_target, type="Cov")
sum(RL2_on_target_cov >= 10)
round(colMeans(RL2_on_target_cov), 1)
RL2_on_target_keep_loci = which(RL2_on_target_cov >= 10)
RL2_on_target_10.X = RL2_on_target[RL2_on_target_keep_loci,]
round(colMeans(getCoverage(RL2_on_target_10.X, type="Cov")), 1)
range(getMeth(RL2_on_target_10.X, type="raw"))
rm(RL2_on_target)
rm(RL2_on_target_10.X)
rm(RL2_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL3.rda")
RL3
RL3_cov = getCoverage(RL3, type="Cov")
sum(RL3_cov >= 10)
round(colMeans(RL3_cov), 1)
RL3_keep_loci = which(RL3_cov >= 10)
RL3_10.X = RL3[RL3_keep_loci,]
round(colMeans(getCoverage(RL3_10.X, type="Cov")), 1)
range(getMeth(RL3_10.X, type="raw"))
rm(RL3)
rm(RL3_10.X)
rm(RL3_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL3_on_target.rda")
RL3_on_target
RL3_on_target_cov = getCoverage(RL3_on_target, type="Cov")
sum(RL3_on_target_cov >= 10)
round(colMeans(RL3_on_target_cov), 1)
RL3_on_target_keep_loci = which(RL3_on_target_cov >= 10)
RL3_on_target_10.X = RL3_on_target[RL3_on_target_keep_loci,]
round(colMeans(getCoverage(RL3_on_target_10.X, type="Cov")), 1)
range(getMeth(RL3_on_target_10.X, type="raw"))
rm(RL3_on_target)
rm(RL3_on_target_10.X)
rm(RL3_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL4.rda")
RL4
RL4_cov = getCoverage(RL4, type="Cov")
sum(RL4_cov >= 10)
round(colMeans(RL4_cov), 1)
RL4_keep_loci = which(RL4_cov >= 10)
RL4_10.X = RL4[RL4_keep_loci,]
round(colMeans(getCoverage(RL4_10.X, type="Cov")), 1)
range(getMeth(RL4_10.X, type="raw"))
rm(RL4)
rm(RL4_10.X)
rm(RL4_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL4_on_target.rda")
RL4_on_target
RL4_on_target_cov = getCoverage(RL4_on_target, type="Cov")
sum(RL4_on_target_cov >= 10)
round(colMeans(RL4_on_target_cov), 1)
RL4_on_target_keep_loci = which(RL4_on_target_cov >= 10)
RL4_on_target_10.X = RL4_on_target[RL4_on_target_keep_loci,]
round(colMeans(getCoverage(RL4_on_target_10.X, type="Cov")), 1)
range(getMeth(RL4_on_target_10.X, type="raw"))
rm(RL4_on_target)
rm(RL4_on_target_10.X)
rm(RL4_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL5.rda")
RL5
RL5_cov = getCoverage(RL5, type="Cov")
sum(RL5_cov >= 10)
round(colMeans(RL5_cov), 1)
RL5_keep_loci = which(RL5_cov >= 10)
RL5_10.X = RL5[RL5_keep_loci,]
round(colMeans(getCoverage(RL5_10.X, type="Cov")), 1)
range(getMeth(RL5_10.X, type="raw"))
rm(RL5)
rm(RL5_10.X)
rm(RL5_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL5_on_target.rda")
RL5_on_target
RL5_on_target_cov = getCoverage(RL5_on_target, type="Cov")
sum(RL5_on_target_cov >= 10)
round(colMeans(RL5_on_target_cov), 1)
RL5_on_target_keep_loci = which(RL5_on_target_cov >= 10)
RL5_on_target_10.X = RL5_on_target[RL5_on_target_keep_loci,]
round(colMeans(getCoverage(RL5_on_target_10.X, type="Cov")), 1)
range(getMeth(RL5_on_target_10.X, type="raw"))
rm(RL5_on_target)
rm(RL5_on_target_10.X)
rm(RL5_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL6.rda")
RL6
RL6_cov = getCoverage(RL6, type="Cov")
sum(RL6_cov >= 10)
round(colMeans(RL6_cov), 1)
RL6_keep_loci = which(RL6_cov >= 10)
RL6_10.X = RL6[RL6_keep_loci,]
round(colMeans(getCoverage(RL6_10.X, type="Cov")), 1)
range(getMeth(RL6_10.X, type="raw"))
rm(RL6)
rm(RL6_10.X)
rm(RL6_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL6_on_target.rda")
RL6_on_target
RL6_on_target_cov = getCoverage(RL6_on_target, type="Cov")
sum(RL6_on_target_cov >= 10)
round(colMeans(RL6_on_target_cov), 1)
RL6_on_target_keep_loci = which(RL6_on_target_cov >= 10)
RL6_on_target_10.X = RL6_on_target[RL6_on_target_keep_loci,]
round(colMeans(getCoverage(RL6_on_target_10.X, type="Cov")), 1)
range(getMeth(RL6_on_target_10.X, type="raw"))
rm(RL6_on_target)
rm(RL6_on_target_10.X)
rm(RL6_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL7.rda")
RL7
RL7_cov = getCoverage(RL7, type="Cov")
sum(RL7_cov >= 10)
round(colMeans(RL7_cov), 1)
RL7_keep_loci = which(RL7_cov >= 10)
RL7_10.X = RL7[RL7_keep_loci,]
round(colMeans(getCoverage(RL7_10.X, type="Cov")), 1)
range(getMeth(RL7_10.X, type="raw"))
rm(RL7)
rm(RL7_10.X)
rm(RL7_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL7_on_target.rda")
RL7_on_target
RL7_on_target_cov = getCoverage(RL7_on_target, type="Cov")
sum(RL7_on_target_cov >= 10)
round(colMeans(RL7_on_target_cov), 1)
RL7_on_target_keep_loci = which(RL7_on_target_cov >= 10)
RL7_on_target_10.X = RL7_on_target[RL7_on_target_keep_loci,]
round(colMeans(getCoverage(RL7_on_target_10.X, type="Cov")), 1)
range(getMeth(RL7_on_target_10.X, type="raw"))
rm(RL7_on_target)
rm(RL7_on_target_10.X)
rm(RL7_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL8.rda")
RL8
RL8_cov = getCoverage(RL8, type="Cov")
sum(RL8_cov >= 10)
round(colMeans(RL8_cov), 1)
RL8_keep_loci = which(RL8_cov >= 10)
RL8_10.X = RL8[RL8_keep_loci,]
round(colMeans(getCoverage(RL8_10.X, type="Cov")), 1)
range(getMeth(RL8_10.X, type="raw"))
rm(RL8)
rm(RL8_10.X)
rm(RL8_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL8_on_target.rda")
RL8_on_target
RL8_on_target_cov = getCoverage(RL8_on_target, type="Cov")
sum(RL8_on_target_cov >= 10)
round(colMeans(RL8_on_target_cov), 1)
RL8_on_target_keep_loci = which(RL8_on_target_cov >= 10)
RL8_on_target_10.X = RL8_on_target[RL8_on_target_keep_loci,]
round(colMeans(getCoverage(RL8_on_target_10.X, type="Cov")), 1)
range(getMeth(RL8_on_target_10.X, type="raw"))
rm(RL8_on_target)
rm(RL8_on_target_10.X)
rm(RL8_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL9.rda")
RL9
RL9_cov = getCoverage(RL9, type="Cov")
sum(RL9_cov >= 10)
round(colMeans(RL9_cov), 1)
RL9_keep_loci = which(RL9_cov >= 10)
RL9_10.X = RL9[RL9_keep_loci,]
round(colMeans(getCoverage(RL9_10.X, type="Cov")), 1)
range(getMeth(RL9_10.X, type="raw"))
rm(RL9)
rm(RL9_10.X)
rm(RL9_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL9_on_target.rda")
RL9_on_target
RL9_on_target_cov = getCoverage(RL9_on_target, type="Cov")
sum(RL9_on_target_cov >= 10)
round(colMeans(RL9_on_target_cov), 1)
RL9_on_target_keep_loci = which(RL9_on_target_cov >= 10)
RL9_on_target_10.X = RL9_on_target[RL9_on_target_keep_loci,]
round(colMeans(getCoverage(RL9_on_target_10.X, type="Cov")), 1)
range(getMeth(RL9_on_target_10.X, type="raw"))
rm(RL9_on_target)
rm(RL9_on_target_10.X)
rm(RL9_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL10.rda")
RL10
RL10_cov = getCoverage(RL10, type="Cov")
sum(RL10_cov >= 10)
round(colMeans(RL10_cov), 1)
RL10_keep_loci = which(RL10_cov >= 10)
RL10_10.X = RL10[RL10_keep_loci,]
round(colMeans(getCoverage(RL10_10.X, type="Cov")), 1)
range(getMeth(RL10_10.X, type="raw"))
rm(RL10)
rm(RL10_10.X)
rm(RL10_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL10_on_target.rda")
RL10_on_target
RL10_on_target_cov = getCoverage(RL10_on_target, type="Cov")
sum(RL10_on_target_cov >= 10)
round(colMeans(RL10_on_target_cov), 1)
RL10_on_target_keep_loci = which(RL10_on_target_cov >= 10)
RL10_on_target_10.X = RL10_on_target[RL10_on_target_keep_loci,]
round(colMeans(getCoverage(RL10_on_target_10.X, type="Cov")), 1)
range(getMeth(RL10_on_target_10.X, type="raw"))
rm(RL10_on_target)
rm(RL10_on_target_10.X)
rm(RL10_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL11.rda")
RL11
RL11_cov = getCoverage(RL11, type="Cov")
sum(RL11_cov >= 10)
round(colMeans(RL11_cov), 1)
RL11_keep_loci = which(RL11_cov >= 10)
RL11_10.X = RL11[RL11_keep_loci,]
round(colMeans(getCoverage(RL11_10.X, type="Cov")), 1)
range(getMeth(RL11_10.X, type="raw"))
rm(RL11)
rm(RL11_10.X)
rm(RL11_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL11_on_target.rda")
RL11_on_target
RL11_on_target_cov = getCoverage(RL11_on_target, type="Cov")
sum(RL11_on_target_cov >= 10)
round(colMeans(RL11_on_target_cov), 1)
RL11_on_target_keep_loci = which(RL11_on_target_cov >= 10)
RL11_on_target_10.X = RL11_on_target[RL11_on_target_keep_loci,]
round(colMeans(getCoverage(RL11_on_target_10.X, type="Cov")), 1)
range(getMeth(RL11_on_target_10.X, type="raw"))
rm(RL11_on_target)
rm(RL11_on_target_10.X)
rm(RL11_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL12.rda")
RL12
RL12_cov = getCoverage(RL12, type="Cov")
sum(RL12_cov >= 10)
round(colMeans(RL12_cov), 1)
RL12_keep_loci = which(RL12_cov >= 10)
RL12_10.X = RL12[RL12_keep_loci,]
round(colMeans(getCoverage(RL12_10.X, type="Cov")), 1)
range(getMeth(RL12_10.X, type="raw"))
rm(RL12)
rm(RL12_10.X)
rm(RL12_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL12_on_target.rda")
RL12_on_target
RL12_on_target_cov = getCoverage(RL12_on_target, type="Cov")
sum(RL12_on_target_cov >= 10)
round(colMeans(RL12_on_target_cov), 1)
RL12_on_target_keep_loci = which(RL12_on_target_cov >= 10)
RL12_on_target_10.X = RL12_on_target[RL12_on_target_keep_loci,]
round(colMeans(getCoverage(RL12_on_target_10.X, type="Cov")), 1)
range(getMeth(RL12_on_target_10.X, type="raw"))
rm(RL12_on_target)
rm(RL12_on_target_10.X)
rm(RL12_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL13.rda")
RL13
RL13_cov = getCoverage(RL13, type="Cov")
sum(RL13_cov >= 10)
round(colMeans(RL13_cov), 1)
RL13_keep_loci = which(RL13_cov >= 10)
RL13_10.X = RL13[RL13_keep_loci,]
round(colMeans(getCoverage(RL13_10.X, type="Cov")), 1)
range(getMeth(RL13_10.X, type="raw"))
rm(RL13)
rm(RL13_10.X)
rm(RL13_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL13_on_target.rda")
RL13_on_target
RL13_on_target_cov = getCoverage(RL13_on_target, type="Cov")
sum(RL13_on_target_cov >= 10)
round(colMeans(RL13_on_target_cov), 1)
RL13_on_target_keep_loci = which(RL13_on_target_cov >= 10)
RL13_on_target_10.X = RL13_on_target[RL13_on_target_keep_loci,]
round(colMeans(getCoverage(RL13_on_target_10.X, type="Cov")), 1)
range(getMeth(RL13_on_target_10.X, type="raw"))
rm(RL13_on_target)
rm(RL13_on_target_10.X)
rm(RL13_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL14.rda")
RL14
RL14_cov = getCoverage(RL14, type="Cov")
sum(RL14_cov >= 10)
round(colMeans(RL14_cov), 1)
RL14_keep_loci = which(RL14_cov >= 10)
RL14_10.X = RL14[RL14_keep_loci,]
round(colMeans(getCoverage(RL14_10.X, type="Cov")), 1)
range(getMeth(RL14_10.X, type="raw"))
rm(RL14)
rm(RL14_10.X)
rm(RL14_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL14_on_target.rda")
RL14_on_target
RL14_on_target_cov = getCoverage(RL14_on_target, type="Cov")
sum(RL14_on_target_cov >= 10)
round(colMeans(RL14_on_target_cov), 1)
RL14_on_target_keep_loci = which(RL14_on_target_cov >= 10)
RL14_on_target_10.X = RL14_on_target[RL14_on_target_keep_loci,]
round(colMeans(getCoverage(RL14_on_target_10.X, type="Cov")), 1)
range(getMeth(RL14_on_target_10.X, type="raw"))
rm(RL14_on_target)
rm(RL14_on_target_10.X)
rm(RL14_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL15.rda")
RL15
RL15_cov = getCoverage(RL15, type="Cov")
sum(RL15_cov >= 10)
round(colMeans(RL15_cov), 1)
RL15_keep_loci = which(RL15_cov >= 10)
RL15_10.X = RL15[RL15_keep_loci,]
round(colMeans(getCoverage(RL15_10.X, type="Cov")), 1)
range(getMeth(RL15_10.X, type="raw"))
rm(RL15)
rm(RL15_10.X)
rm(RL15_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL15_on_target.rda")
RL15_on_target
RL15_on_target_cov = getCoverage(RL15_on_target, type="Cov")
sum(RL15_on_target_cov >= 10)
round(colMeans(RL15_on_target_cov), 1)
RL15_on_target_keep_loci = which(RL15_on_target_cov >= 10)
RL15_on_target_10.X = RL15_on_target[RL15_on_target_keep_loci,]
round(colMeans(getCoverage(RL15_on_target_10.X, type="Cov")), 1)
range(getMeth(RL15_on_target_10.X, type="raw"))
rm(RL15_on_target)
rm(RL15_on_target_10.X)
rm(RL15_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL16.rda")
RL16
RL16_cov = getCoverage(RL16, type="Cov")
sum(RL16_cov >= 10)
round(colMeans(RL16_cov), 1)
RL16_keep_loci = which(RL16_cov >= 10)
RL16_10.X = RL16[RL16_keep_loci,]
round(colMeans(getCoverage(RL16_10.X, type="Cov")), 1)
range(getMeth(RL16_10.X, type="raw"))
rm(RL16)
rm(RL16_10.X)
rm(RL16_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL16_on_target.rda")
RL16_on_target
RL16_on_target_cov = getCoverage(RL16_on_target, type="Cov")
sum(RL16_on_target_cov >= 10)
round(colMeans(RL16_on_target_cov), 1)
RL16_on_target_keep_loci = which(RL16_on_target_cov >= 10)
RL16_on_target_10.X = RL16_on_target[RL16_on_target_keep_loci,]
round(colMeans(getCoverage(RL16_on_target_10.X, type="Cov")), 1)
range(getMeth(RL16_on_target_10.X, type="raw"))
rm(RL16_on_target)
rm(RL16_on_target_10.X)
rm(RL16_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL17.rda")
RL17
RL17_cov = getCoverage(RL17, type="Cov")
sum(RL17_cov >= 10)
round(colMeans(RL17_cov), 1)
RL17_keep_loci = which(RL17_cov >= 10)
RL17_10.X = RL17[RL17_keep_loci,]
round(colMeans(getCoverage(RL17_10.X, type="Cov")), 1)
range(getMeth(RL17_10.X, type="raw"))
rm(RL17)
rm(RL17_10.X)
rm(RL17_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL17_on_target.rda")
RL17_on_target
RL17_on_target_cov = getCoverage(RL17_on_target, type="Cov")
sum(RL17_on_target_cov >= 10)
round(colMeans(RL17_on_target_cov), 1)
RL17_on_target_keep_loci = which(RL17_on_target_cov >= 10)
RL17_on_target_10.X = RL17_on_target[RL17_on_target_keep_loci,]
round(colMeans(getCoverage(RL17_on_target_10.X, type="Cov")), 1)
range(getMeth(RL17_on_target_10.X, type="raw"))
rm(RL17_on_target)
rm(RL17_on_target_10.X)
rm(RL17_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL18.rda")
RL18
RL18_cov = getCoverage(RL18, type="Cov")
sum(RL18_cov >= 10)
round(colMeans(RL18_cov), 1)
RL18_keep_loci = which(RL18_cov >= 10)
RL18_10.X = RL18[RL18_keep_loci,]
round(colMeans(getCoverage(RL18_10.X, type="Cov")), 1)
range(getMeth(RL18_10.X, type="raw"))
rm(RL18)
rm(RL18_10.X)
rm(RL18_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL18_on_target.rda")
RL18_on_target
RL18_on_target_cov = getCoverage(RL18_on_target, type="Cov")
sum(RL18_on_target_cov >= 10)
round(colMeans(RL18_on_target_cov), 1)
RL18_on_target_keep_loci = which(RL18_on_target_cov >= 10)
RL18_on_target_10.X = RL18_on_target[RL18_on_target_keep_loci,]
round(colMeans(getCoverage(RL18_on_target_10.X, type="Cov")), 1)
range(getMeth(RL18_on_target_10.X, type="raw"))
rm(RL18_on_target)
rm(RL18_on_target_10.X)
rm(RL18_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL19.rda")
RL19
RL19_cov = getCoverage(RL19, type="Cov")
sum(RL19_cov >= 10)
round(colMeans(RL19_cov), 1)
RL19_keep_loci = which(RL19_cov >= 10)
RL19_10.X = RL19[RL19_keep_loci,]
round(colMeans(getCoverage(RL19_10.X, type="Cov")), 1)
range(getMeth(RL19_10.X, type="raw"))
rm(RL19)
rm(RL19_10.X)
rm(RL19_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL19_on_target.rda")
RL19_on_target
RL19_on_target_cov = getCoverage(RL19_on_target, type="Cov")
sum(RL19_on_target_cov >= 10)
round(colMeans(RL19_on_target_cov), 1)
RL19_on_target_keep_loci = which(RL19_on_target_cov >= 10)
RL19_on_target_10.X = RL19_on_target[RL19_on_target_keep_loci,]
round(colMeans(getCoverage(RL19_on_target_10.X, type="Cov")), 1)
range(getMeth(RL19_on_target_10.X, type="raw"))
rm(RL19_on_target)
rm(RL19_on_target_10.X)
rm(RL19_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL20.rda")
RL20
RL20_cov = getCoverage(RL20, type="Cov")
sum(RL20_cov >= 10)
round(colMeans(RL20_cov), 1)
RL20_keep_loci = which(RL20_cov >= 10)
RL20_10.X = RL20[RL20_keep_loci,]
round(colMeans(getCoverage(RL20_10.X, type="Cov")), 1)
range(getMeth(RL20_10.X, type="raw"))
rm(RL20)
rm(RL20_10.X)
rm(RL20_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL20_on_target.rda")
RL20_on_target
RL20_on_target_cov = getCoverage(RL20_on_target, type="Cov")
sum(RL20_on_target_cov >= 10)
round(colMeans(RL20_on_target_cov), 1)
RL20_on_target_keep_loci = which(RL20_on_target_cov >= 10)
RL20_on_target_10.X = RL20_on_target[RL20_on_target_keep_loci,]
round(colMeans(getCoverage(RL20_on_target_10.X, type="Cov")), 1)
range(getMeth(RL20_on_target_10.X, type="raw"))
rm(RL20_on_target)
rm(RL20_on_target_10.X)
rm(RL20_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL21.rda")
RL21
RL21_cov = getCoverage(RL21, type="Cov")
sum(RL21_cov >= 10)
round(colMeans(RL21_cov), 1)
RL21_keep_loci = which(RL21_cov >= 10)
RL21_10.X = RL21[RL21_keep_loci,]
round(colMeans(getCoverage(RL21_10.X, type="Cov")), 1)
range(getMeth(RL21_10.X, type="raw"))
rm(RL21)
rm(RL21_10.X)
rm(RL21_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL21_on_target.rda")
RL21_on_target
RL21_on_target_cov = getCoverage(RL21_on_target, type="Cov")
sum(RL21_on_target_cov >= 10)
round(colMeans(RL21_on_target_cov), 1)
RL21_on_target_keep_loci = which(RL21_on_target_cov >= 10)
RL21_on_target_10.X = RL21_on_target[RL21_on_target_keep_loci,]
round(colMeans(getCoverage(RL21_on_target_10.X, type="Cov")), 1)
range(getMeth(RL21_on_target_10.X, type="raw"))
rm(RL21_on_target)
rm(RL21_on_target_10.X)
rm(RL21_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL22.rda")
RL22
RL22_cov = getCoverage(RL22, type="Cov")
sum(RL22_cov >= 10)
round(colMeans(RL22_cov), 1)
RL22_keep_loci = which(RL22_cov >= 10)
RL22_10.X = RL22[RL22_keep_loci,]
round(colMeans(getCoverage(RL22_10.X, type="Cov")), 1)
range(getMeth(RL22_10.X, type="raw"))
rm(RL22)
rm(RL22_10.X)
rm(RL22_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL22_on_target.rda")
RL22_on_target
RL22_on_target_cov = getCoverage(RL22_on_target, type="Cov")
sum(RL22_on_target_cov >= 10)
round(colMeans(RL22_on_target_cov), 1)
RL22_on_target_keep_loci = which(RL22_on_target_cov >= 10)
RL22_on_target_10.X = RL22_on_target[RL22_on_target_keep_loci,]
round(colMeans(getCoverage(RL22_on_target_10.X, type="Cov")), 1)
range(getMeth(RL22_on_target_10.X, type="raw"))
rm(RL22_on_target)
rm(RL22_on_target_10.X)
rm(RL22_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL23.rda")
RL23
RL23_cov = getCoverage(RL23, type="Cov")
sum(RL23_cov >= 10)
round(colMeans(RL23_cov), 1)
RL23_keep_loci = which(RL23_cov >= 10)
RL23_10.X = RL23[RL23_keep_loci,]
round(colMeans(getCoverage(RL23_10.X, type="Cov")), 1)
range(getMeth(RL23_10.X, type="raw"))
rm(RL23)
rm(RL23_10.X)
rm(RL23_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL23_on_target.rda")
RL23_on_target
RL23_on_target_cov = getCoverage(RL23_on_target, type="Cov")
sum(RL23_on_target_cov >= 10)
round(colMeans(RL23_on_target_cov), 1)
RL23_on_target_keep_loci = which(RL23_on_target_cov >= 10)
RL23_on_target_10.X = RL23_on_target[RL23_on_target_keep_loci,]
round(colMeans(getCoverage(RL23_on_target_10.X, type="Cov")), 1)
range(getMeth(RL23_on_target_10.X, type="raw"))
rm(RL23_on_target)
rm(RL23_on_target_10.X)
rm(RL23_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL24.rda")
RL24
RL24_cov = getCoverage(RL24, type="Cov")
sum(RL24_cov >= 10)
round(colMeans(RL24_cov), 1)
RL24_keep_loci = which(RL24_cov >= 10)
RL24_10.X = RL24[RL24_keep_loci,]
round(colMeans(getCoverage(RL24_10.X, type="Cov")), 1)
range(getMeth(RL24_10.X, type="raw"))
rm(RL24)
rm(RL24_10.X)
rm(RL24_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL24_on_target.rda")
RL24_on_target
RL24_on_target_cov = getCoverage(RL24_on_target, type="Cov")
sum(RL24_on_target_cov >= 10)
round(colMeans(RL24_on_target_cov), 1)
RL24_on_target_keep_loci = which(RL24_on_target_cov >= 10)
RL24_on_target_10.X = RL24_on_target[RL24_on_target_keep_loci,]
round(colMeans(getCoverage(RL24_on_target_10.X, type="Cov")), 1)
range(getMeth(RL24_on_target_10.X, type="raw"))
rm(RL24_on_target)
rm(RL24_on_target_10.X)
rm(RL24_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL25.rda")
RL25
RL25_cov = getCoverage(RL25, type="Cov")
sum(RL25_cov >= 10)
round(colMeans(RL25_cov), 1)
RL25_keep_loci = which(RL25_cov >= 10)
RL25_10.X = RL25[RL25_keep_loci,]
round(colMeans(getCoverage(RL25_10.X, type="Cov")), 1)
range(getMeth(RL25_10.X, type="raw"))
rm(RL25)
rm(RL25_10.X)
rm(RL25_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL25_on_target.rda")
RL25_on_target
RL25_on_target_cov = getCoverage(RL25_on_target, type="Cov")
sum(RL25_on_target_cov >= 10)
round(colMeans(RL25_on_target_cov), 1)
RL25_on_target_keep_loci = which(RL25_on_target_cov >= 10)
RL25_on_target_10.X = RL25_on_target[RL25_on_target_keep_loci,]
round(colMeans(getCoverage(RL25_on_target_10.X, type="Cov")), 1)
range(getMeth(RL25_on_target_10.X, type="raw"))
rm(RL25_on_target)
rm(RL25_on_target_10.X)
rm(RL25_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL26.rda")
RL26
RL26_cov = getCoverage(RL26, type="Cov")
sum(RL26_cov >= 10)
round(colMeans(RL26_cov), 1)
RL26_keep_loci = which(RL26_cov >= 10)
RL26_10.X = RL26[RL26_keep_loci,]
round(colMeans(getCoverage(RL26_10.X, type="Cov")), 1)
range(getMeth(RL26_10.X, type="raw"))
rm(RL26)
rm(RL26_10.X)
rm(RL26_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL26_on_target.rda")
RL26_on_target
RL26_on_target_cov = getCoverage(RL26_on_target, type="Cov")
sum(RL26_on_target_cov >= 10)
round(colMeans(RL26_on_target_cov), 1)
RL26_on_target_keep_loci = which(RL26_on_target_cov >= 10)
RL26_on_target_10.X = RL26_on_target[RL26_on_target_keep_loci,]
round(colMeans(getCoverage(RL26_on_target_10.X, type="Cov")), 1)
range(getMeth(RL26_on_target_10.X, type="raw"))
rm(RL26_on_target)
rm(RL26_on_target_10.X)
rm(RL26_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL27.rda")
RL27
RL27_cov = getCoverage(RL27, type="Cov")
sum(RL27_cov >= 10)
round(colMeans(RL27_cov), 1)
RL27_keep_loci = which(RL27_cov >= 10)
RL27_10.X = RL27[RL27_keep_loci,]
round(colMeans(getCoverage(RL27_10.X, type="Cov")), 1)
range(getMeth(RL27_10.X, type="raw"))
rm(RL27)
rm(RL27_10.X)
rm(RL27_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL27_on_target.rda")
RL27_on_target
RL27_on_target_cov = getCoverage(RL27_on_target, type="Cov")
sum(RL27_on_target_cov >= 10)
round(colMeans(RL27_on_target_cov), 1)
RL27_on_target_keep_loci = which(RL27_on_target_cov >= 10)
RL27_on_target_10.X = RL27_on_target[RL27_on_target_keep_loci,]
round(colMeans(getCoverage(RL27_on_target_10.X, type="Cov")), 1)
range(getMeth(RL27_on_target_10.X, type="raw"))
rm(RL27_on_target)
rm(RL27_on_target_10.X)
rm(RL27_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL28.rda")
RL28
RL28_cov = getCoverage(RL28, type="Cov")
sum(RL28_cov >= 10)
round(colMeans(RL28_cov), 1)
RL28_keep_loci = which(RL28_cov >= 10)
RL28_10.X = RL28[RL28_keep_loci,]
round(colMeans(getCoverage(RL28_10.X, type="Cov")), 1)
range(getMeth(RL28_10.X, type="raw"))
rm(RL28)
rm(RL28_10.X)
rm(RL28_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL28_on_target.rda")
RL28_on_target
RL28_on_target_cov = getCoverage(RL28_on_target, type="Cov")
sum(RL28_on_target_cov >= 10)
round(colMeans(RL28_on_target_cov), 1)
RL28_on_target_keep_loci = which(RL28_on_target_cov >= 10)
RL28_on_target_10.X = RL28_on_target[RL28_on_target_keep_loci,]
round(colMeans(getCoverage(RL28_on_target_10.X, type="Cov")), 1)
range(getMeth(RL28_on_target_10.X, type="raw"))
rm(RL28_on_target)
rm(RL28_on_target_10.X)
rm(RL28_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL29.rda")
RL29
RL29_cov = getCoverage(RL29, type="Cov")
sum(RL29_cov >= 10)
round(colMeans(RL29_cov), 1)
RL29_keep_loci = which(RL29_cov >= 10)
RL29_10.X = RL29[RL29_keep_loci,]
round(colMeans(getCoverage(RL29_10.X, type="Cov")), 1)
range(getMeth(RL29_10.X, type="raw"))
rm(RL29)
rm(RL29_10.X)
rm(RL29_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL29_on_target.rda")
RL29_on_target
RL29_on_target_cov = getCoverage(RL29_on_target, type="Cov")
sum(RL29_on_target_cov >= 10)
round(colMeans(RL29_on_target_cov), 1)
RL29_on_target_keep_loci = which(RL29_on_target_cov >= 10)
RL29_on_target_10.X = RL29_on_target[RL29_on_target_keep_loci,]
round(colMeans(getCoverage(RL29_on_target_10.X, type="Cov")), 1)
range(getMeth(RL29_on_target_10.X, type="raw"))
rm(RL29_on_target)
rm(RL29_on_target_10.X)
rm(RL29_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL30.rda")
RL30
RL30_cov = getCoverage(RL30, type="Cov")
sum(RL30_cov >= 10)
round(colMeans(RL30_cov), 1)
RL30_keep_loci = which(RL30_cov >= 10)
RL30_10.X = RL30[RL30_keep_loci,]
round(colMeans(getCoverage(RL30_10.X, type="Cov")), 1)
range(getMeth(RL30_10.X, type="raw"))
rm(RL30)
rm(RL30_10.X)
rm(RL30_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL30_on_target.rda")
RL30_on_target
RL30_on_target_cov = getCoverage(RL30_on_target, type="Cov")
sum(RL30_on_target_cov >= 10)
round(colMeans(RL30_on_target_cov), 1)
RL30_on_target_keep_loci = which(RL30_on_target_cov >= 10)
RL30_on_target_10.X = RL30_on_target[RL30_on_target_keep_loci,]
round(colMeans(getCoverage(RL30_on_target_10.X, type="Cov")), 1)
range(getMeth(RL30_on_target_10.X, type="raw"))
rm(RL30_on_target)
rm(RL30_on_target_10.X)
rm(RL30_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL31.rda")
RL31
RL31_cov = getCoverage(RL31, type="Cov")
sum(RL31_cov >= 10)
round(colMeans(RL31_cov), 1)
RL31_keep_loci = which(RL31_cov >= 10)
RL31_10.X = RL31[RL31_keep_loci,]
round(colMeans(getCoverage(RL31_10.X, type="Cov")), 1)
range(getMeth(RL31_10.X, type="raw"))
rm(RL31)
rm(RL31_10.X)
rm(RL31_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL31_on_target.rda")
RL31_on_target
RL31_on_target_cov = getCoverage(RL31_on_target, type="Cov")
sum(RL31_on_target_cov >= 10)
round(colMeans(RL31_on_target_cov), 1)
RL31_on_target_keep_loci = which(RL31_on_target_cov >= 10)
RL31_on_target_10.X = RL31_on_target[RL31_on_target_keep_loci,]
round(colMeans(getCoverage(RL31_on_target_10.X, type="Cov")), 1)
range(getMeth(RL31_on_target_10.X, type="raw"))
rm(RL31_on_target)
rm(RL31_on_target_10.X)
rm(RL31_on_target_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL32.rda")
RL32
RL32_cov = getCoverage(RL32, type="Cov")
sum(RL32_cov >= 10)
round(colMeans(RL32_cov), 1)
RL32_keep_loci = which(RL32_cov >= 10)
RL32_10.X = RL32[RL32_keep_loci,]
round(colMeans(getCoverage(RL32_10.X, type="Cov")), 1)
range(getMeth(RL32_10.X, type="raw"))
rm(RL32)
rm(RL32_10.X)
rm(RL32_cov)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL32_on_target.rda")
RL32_on_target
RL32_on_target_cov = getCoverage(RL32_on_target, type="Cov")
sum(RL32_on_target_cov >= 10)
round(colMeans(RL32_on_target_cov), 1)
RL32_on_target_keep_loci = which(RL32_on_target_cov >= 10)
RL32_on_target_10.X = RL32_on_target[RL32_on_target_keep_loci,]
round(colMeans(getCoverage(RL32_on_target_10.X, type="Cov")), 1)
range(getMeth(RL32_on_target_10.X, type="raw"))
rm(RL32_on_target)
rm(RL32_on_target_10.X)
rm(RL32_on_target_cov)

