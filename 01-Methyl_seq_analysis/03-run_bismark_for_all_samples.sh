#!/bin/bash
#$ -cwd

datadir="/home/cluster/nas1/fseifudd/02-methylseq/data-source/gary_stress_human_methylseq"
workingdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/04-alignment_with_bismark"

while read sampleid
do
	mkdir $workingdir/$sampleid
	sh /home/cluster/nas1/fseifudd/02-methylseq/methylseq-box/run_bismark.sh $datadir/JHURL01009_$sampleid/"$sampleid"_R1.fastq.gz $datadir/JHURL01009_$sampleid/"$sampleid"_R2.fastq.gz $sampleid
done < $workingdir/gary_human_methylseq_sample_ids.txt
