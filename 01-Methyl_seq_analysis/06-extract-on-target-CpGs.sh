#!/bin/bash
#$ -cwd

. /home/cluster/nas1/fseifudd/02-methylseq/methylseq-box/config.sh

datadir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark"
outputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/08-generate-on-target-CpGs"
bedfile="/home/cluster/nas1/fseifudd/02-methylseq/data-source"

while read subjectid
do
	bedtools intersect -a $bedfile/SureSelect_Human_Methyl_Seq_Kit.bed -b $datadir/$subjectid/bismark_out/methylation_extractor_output_deduplicated/$subjectid.bismark_pe.deduplicated.bismark.cov.gz > $datadir/$subjectid/bismark_out/methylation_extractor_output_deduplicated/$subjectid.bismark_pe.deduplicated.bismark.on.target.cov
	
done < $outputdir/gary_human_methylseq_sample_ids.txt

