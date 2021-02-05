#!/bin/bash
#$ -cwd

. /home/cluster/nas1/fseifudd/02-methylseq/methylseq-box/config.sh

datadir="/home/cluster/nas1/fseifudd/02-methylseq/data-source/gary_stress_human_methylseq"
outputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/01-fastqc"

while read sampleid
do
	echo processing $sampleid ...
#	mkdir $outputdir/$sampleid
	perl $fastqc_dir/fastqc -o $outputdir/$sampleid $datadir/JHURL01009_$sampleid/"$sampleid"_R1.fastq.gz $datadir/JHURL01009_$sampleid/"$sampleid"_R2.fastq.gz
	echo " "
	echo " "
	echo finished proccessing $sampleid
done < $outputdir/gary_human_methylseq_sample_ids.txt

