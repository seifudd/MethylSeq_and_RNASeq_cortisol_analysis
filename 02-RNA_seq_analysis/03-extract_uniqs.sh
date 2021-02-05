#!/bin/bash

# exclude unmapped (-F 4) 
# exclude multiple-mapped (grep -v "XS:") 

bam_dir="$1"
bam_f="$2"

module load samtools

samtools view -H   $bam_dir/$bam_f  > $bam_dir/$bam_f.header.sam
samtools view -f 3 $bam_dir/$bam_f | grep "NH:i:1"  | cat $bam_dir/$bam_f.header.sam -  \
 | samtools view -b - > $bam_dir/$bam_f.unique.bam

samtools index $bam_dir/$bam_f.unique.bam

rm $bam_dir/$bam_f.header.sam


