#!/bin/bash
#$ -cwd

inputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark"
workingdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/06-methylation_extraction_CpG_context"

while  read subjectid
do
#	mkdir $workingdir/$subjectid
	sh /home/cluster/nas1/fseifudd/02-methylseq/methylseq-box/run_bismark_methylation_extractor.sh $subjectid $inputdir $workingdir
done < $workingdir/gary_human_methylseq_sample_ids.txt

