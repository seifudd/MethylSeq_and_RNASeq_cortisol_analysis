library(bsseq)

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL1.rda")
RL1_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL1/bismark_out/methylation_extractor_output_deduplicated/RL1.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL1_on_target = subsetByOverlaps(RL1, GRanges(seqnames=RL1_on_target_cpgs[,1], ranges= IRanges(start=RL1_on_target_cpgs[,2], end=RL1_on_target_cpgs[,3])))
RL1_on_target
save(RL1_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL1_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL2.rda")
RL2_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL2/bismark_out/methylation_extractor_output_deduplicated/RL2.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL2_on_target = subsetByOverlaps(RL2, GRanges(seqnames=RL2_on_target_cpgs[,1], ranges= IRanges(start=RL2_on_target_cpgs[,2], end=RL2_on_target_cpgs[,3])))
RL2_on_target
save(RL2_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL2_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL3.rda")
RL3_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL3/bismark_out/methylation_extractor_output_deduplicated/RL3.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL3_on_target = subsetByOverlaps(RL3, GRanges(seqnames=RL3_on_target_cpgs[,1], ranges= IRanges(start=RL3_on_target_cpgs[,2], end=RL3_on_target_cpgs[,3])))
RL3_on_target
save(RL3_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL3_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL4.rda")
RL4_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL4/bismark_out/methylation_extractor_output_deduplicated/RL4.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL4_on_target = subsetByOverlaps(RL4, GRanges(seqnames=RL4_on_target_cpgs[,1], ranges= IRanges(start=RL4_on_target_cpgs[,2], end=RL4_on_target_cpgs[,3])))
RL4_on_target
save(RL4_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL4_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL5.rda")
RL5_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL5/bismark_out/methylation_extractor_output_deduplicated/RL5.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL5_on_target = subsetByOverlaps(RL5, GRanges(seqnames=RL5_on_target_cpgs[,1], ranges= IRanges(start=RL5_on_target_cpgs[,2], end=RL5_on_target_cpgs[,3])))
RL5_on_target
save(RL5_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL5_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL6.rda")
RL6_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL6/bismark_out/methylation_extractor_output_deduplicated/RL6.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL6_on_target = subsetByOverlaps(RL6, GRanges(seqnames=RL6_on_target_cpgs[,1], ranges= IRanges(start=RL6_on_target_cpgs[,2], end=RL6_on_target_cpgs[,3])))
RL6_on_target
save(RL6_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL6_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL7.rda")
RL7_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL7/bismark_out/methylation_extractor_output_deduplicated/RL7.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL7_on_target = subsetByOverlaps(RL7, GRanges(seqnames=RL7_on_target_cpgs[,1], ranges= IRanges(start=RL7_on_target_cpgs[,2], end=RL7_on_target_cpgs[,3])))
RL7_on_target
save(RL7_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL7_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL8.rda")
RL8_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL8/bismark_out/methylation_extractor_output_deduplicated/RL8.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL8_on_target = subsetByOverlaps(RL8, GRanges(seqnames=RL8_on_target_cpgs[,1], ranges= IRanges(start=RL8_on_target_cpgs[,2], end=RL8_on_target_cpgs[,3])))
RL8_on_target
save(RL8_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL8_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL9.rda")
RL9_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL9/bismark_out/methylation_extractor_output_deduplicated/RL9.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL9_on_target = subsetByOverlaps(RL9, GRanges(seqnames=RL9_on_target_cpgs[,1], ranges= IRanges(start=RL9_on_target_cpgs[,2], end=RL9_on_target_cpgs[,3])))
RL9_on_target
save(RL9_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL9_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL10.rda")
RL10_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL10/bismark_out/methylation_extractor_output_deduplicated/RL10.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL10_on_target = subsetByOverlaps(RL10, GRanges(seqnames=RL10_on_target_cpgs[,1], ranges= IRanges(start=RL10_on_target_cpgs[,2], end=RL10_on_target_cpgs[,3])))
RL10_on_target
save(RL10_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL10_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL11.rda")
RL11_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL11/bismark_out/methylation_extractor_output_deduplicated/RL11.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL11_on_target = subsetByOverlaps(RL11, GRanges(seqnames=RL11_on_target_cpgs[,1], ranges= IRanges(start=RL11_on_target_cpgs[,2], end=RL11_on_target_cpgs[,3])))
RL11_on_target
save(RL11_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL11_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL12.rda")
RL12_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL12/bismark_out/methylation_extractor_output_deduplicated/RL12.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL12_on_target = subsetByOverlaps(RL12, GRanges(seqnames=RL12_on_target_cpgs[,1], ranges= IRanges(start=RL12_on_target_cpgs[,2], end=RL12_on_target_cpgs[,3])))
RL12_on_target
save(RL12_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL12_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL13.rda")
RL13_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL13/bismark_out/methylation_extractor_output_deduplicated/RL13.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL13_on_target = subsetByOverlaps(RL13, GRanges(seqnames=RL13_on_target_cpgs[,1], ranges= IRanges(start=RL13_on_target_cpgs[,2], end=RL13_on_target_cpgs[,3])))
RL13_on_target
save(RL13_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL13_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL14.rda")
RL14_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL14/bismark_out/methylation_extractor_output_deduplicated/RL14.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL14_on_target = subsetByOverlaps(RL14, GRanges(seqnames=RL14_on_target_cpgs[,1], ranges= IRanges(start=RL14_on_target_cpgs[,2], end=RL14_on_target_cpgs[,3])))
RL14_on_target
save(RL14_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL14_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL15.rda")
RL15_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL15/bismark_out/methylation_extractor_output_deduplicated/RL15.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL15_on_target = subsetByOverlaps(RL15, GRanges(seqnames=RL15_on_target_cpgs[,1], ranges= IRanges(start=RL15_on_target_cpgs[,2], end=RL15_on_target_cpgs[,3])))
RL15_on_target
save(RL15_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL15_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL16.rda")
RL16_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL16/bismark_out/methylation_extractor_output_deduplicated/RL16.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL16_on_target = subsetByOverlaps(RL16, GRanges(seqnames=RL16_on_target_cpgs[,1], ranges= IRanges(start=RL16_on_target_cpgs[,2], end=RL16_on_target_cpgs[,3])))
RL16_on_target
save(RL16_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL16_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL17.rda")
RL17_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL17/bismark_out/methylation_extractor_output_deduplicated/RL17.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL17_on_target = subsetByOverlaps(RL17, GRanges(seqnames=RL17_on_target_cpgs[,1], ranges= IRanges(start=RL17_on_target_cpgs[,2], end=RL17_on_target_cpgs[,3])))
RL17_on_target
save(RL17_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL17_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL18.rda")
RL18_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL18/bismark_out/methylation_extractor_output_deduplicated/RL18.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL18_on_target = subsetByOverlaps(RL18, GRanges(seqnames=RL18_on_target_cpgs[,1], ranges= IRanges(start=RL18_on_target_cpgs[,2], end=RL18_on_target_cpgs[,3])))
RL18_on_target
save(RL18_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL18_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL19.rda")
RL19_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL19/bismark_out/methylation_extractor_output_deduplicated/RL19.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL19_on_target = subsetByOverlaps(RL19, GRanges(seqnames=RL19_on_target_cpgs[,1], ranges= IRanges(start=RL19_on_target_cpgs[,2], end=RL19_on_target_cpgs[,3])))
RL19_on_target
save(RL19_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL19_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL20.rda")
RL20_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL20/bismark_out/methylation_extractor_output_deduplicated/RL20.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL20_on_target = subsetByOverlaps(RL20, GRanges(seqnames=RL20_on_target_cpgs[,1], ranges= IRanges(start=RL20_on_target_cpgs[,2], end=RL20_on_target_cpgs[,3])))
RL20_on_target
save(RL20_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL20_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL21.rda")
RL21_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL21/bismark_out/methylation_extractor_output_deduplicated/RL21.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL21_on_target = subsetByOverlaps(RL21, GRanges(seqnames=RL21_on_target_cpgs[,1], ranges= IRanges(start=RL21_on_target_cpgs[,2], end=RL21_on_target_cpgs[,3])))
RL21_on_target
save(RL21_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL21_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL22.rda")
RL22_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL22/bismark_out/methylation_extractor_output_deduplicated/RL22.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL22_on_target = subsetByOverlaps(RL22, GRanges(seqnames=RL22_on_target_cpgs[,1], ranges= IRanges(start=RL22_on_target_cpgs[,2], end=RL22_on_target_cpgs[,3])))
RL22_on_target
save(RL22_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL22_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL23.rda")
RL23_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL23/bismark_out/methylation_extractor_output_deduplicated/RL23.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL23_on_target = subsetByOverlaps(RL23, GRanges(seqnames=RL23_on_target_cpgs[,1], ranges= IRanges(start=RL23_on_target_cpgs[,2], end=RL23_on_target_cpgs[,3])))
RL23_on_target
save(RL23_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL23_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL24.rda")
RL24_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL24/bismark_out/methylation_extractor_output_deduplicated/RL24.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL24_on_target = subsetByOverlaps(RL24, GRanges(seqnames=RL24_on_target_cpgs[,1], ranges= IRanges(start=RL24_on_target_cpgs[,2], end=RL24_on_target_cpgs[,3])))
RL24_on_target
save(RL24_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL24_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL25.rda")
RL25_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL25/bismark_out/methylation_extractor_output_deduplicated/RL25.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL25_on_target = subsetByOverlaps(RL25, GRanges(seqnames=RL25_on_target_cpgs[,1], ranges= IRanges(start=RL25_on_target_cpgs[,2], end=RL25_on_target_cpgs[,3])))
RL25_on_target
save(RL25_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL25_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL26.rda")
RL26_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL26/bismark_out/methylation_extractor_output_deduplicated/RL26.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL26_on_target = subsetByOverlaps(RL26, GRanges(seqnames=RL26_on_target_cpgs[,1], ranges= IRanges(start=RL26_on_target_cpgs[,2], end=RL26_on_target_cpgs[,3])))
RL26_on_target
save(RL26_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL26_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL27.rda")
RL27_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL27/bismark_out/methylation_extractor_output_deduplicated/RL27.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL27_on_target = subsetByOverlaps(RL27, GRanges(seqnames=RL27_on_target_cpgs[,1], ranges= IRanges(start=RL27_on_target_cpgs[,2], end=RL27_on_target_cpgs[,3])))
RL27_on_target
save(RL27_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL27_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL28.rda")
RL28_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL28/bismark_out/methylation_extractor_output_deduplicated/RL28.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL28_on_target = subsetByOverlaps(RL28, GRanges(seqnames=RL28_on_target_cpgs[,1], ranges= IRanges(start=RL28_on_target_cpgs[,2], end=RL28_on_target_cpgs[,3])))
RL28_on_target
save(RL28_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL28_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL29.rda")
RL29_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL29/bismark_out/methylation_extractor_output_deduplicated/RL29.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL29_on_target = subsetByOverlaps(RL29, GRanges(seqnames=RL29_on_target_cpgs[,1], ranges= IRanges(start=RL29_on_target_cpgs[,2], end=RL29_on_target_cpgs[,3])))
RL29_on_target
save(RL29_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL29_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL30.rda")
RL30_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL30/bismark_out/methylation_extractor_output_deduplicated/RL30.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL30_on_target = subsetByOverlaps(RL30, GRanges(seqnames=RL30_on_target_cpgs[,1], ranges= IRanges(start=RL30_on_target_cpgs[,2], end=RL30_on_target_cpgs[,3])))
RL30_on_target
save(RL30_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL30_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL31.rda")
RL31_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL31/bismark_out/methylation_extractor_output_deduplicated/RL31.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL31_on_target = subsetByOverlaps(RL31, GRanges(seqnames=RL31_on_target_cpgs[,1], ranges= IRanges(start=RL31_on_target_cpgs[,2], end=RL31_on_target_cpgs[,3])))
RL31_on_target
save(RL31_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL31_on_target.rda")

load("/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/09-bsseq_R_data_object_by_sample/RL32.rda")
RL32_on_target_cpgs = read.table(file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark/RL32/bismark_out/methylation_extractor_output_deduplicated/RL32.bismark_pe.deduplicated.bismark.on.target.cov", header=F)
RL32_on_target = subsetByOverlaps(RL32, GRanges(seqnames=RL32_on_target_cpgs[,1], ranges= IRanges(start=RL32_on_target_cpgs[,2], end=RL32_on_target_cpgs[,3])))
RL32_on_target
save(RL32_on_target, file="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/10-bsseq_R_data_object_by_sample_on_target/RL32_on_target.rda")

