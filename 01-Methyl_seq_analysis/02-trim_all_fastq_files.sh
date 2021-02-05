#!/bin/bash
#$ -cwd

. /home/cluster/nas1/fseifudd/02-methylseq/methylseq-box/config.sh

datadir="/home/cluster/nas1/fseifudd/02-methylseq/data-source/gary_stress_human_methylseq"
outputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/02-trim_galore"
fastqc_outputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/03-fastqc_after_trimgalore"

while read sampleid
do
	echo processing $sampleid ...
	mkdir $outputdir/$sampleid
	mkdir $fastqc_outputdir/$sampleid
	perl $trimgalore_dir/trim_galore -o $outputdir/$sampleid --stringency 5 --fastqc_args "-o $fastqc_outputdir/$sampleid" -path_to_cutadapt="/home/cluster/nas1/fseifudd/02-methylseq/application/cutadapt-1.8.1/bin/cutadapt" -paired $datadir/JHURL01009_$sampleid/"$sampleid"_R1.fastq.gz $datadir/JHURL01009_$sampleid/"$sampleid"_R2.fastq.gz 
	echo " "
	echo " "
	echo finished proccessing $sampleid
done < $outputdir/gary_human_methylseq_sample_ids.txt
