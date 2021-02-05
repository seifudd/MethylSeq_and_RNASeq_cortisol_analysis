#!/bin/bash
set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

fq_1="$1"
fq_2="$2"
sid="$3"
pdir="$4"
REF="$5"

output_d="$pdir/03-alignment/$sid"

module load hisat
module load samtools/1.3

mkdir -p $output_d

hisat2 -p 8 \
  -x $REF \
  --dta-cufflinks \
  -1 $fq_1 \
  -2 $fq_2 \
  -S $output_d/$sid.out.sam \
  --rg-id $sid --rg SM:$sid 

samtools view -b -S $output_d/$sid.out.sam > $output_d/$sid.out.bam

samtools sort $output_d/$sid.out.bam -o $output_d/$sid.out.sorted.bam

samtools index $output_d/$sid.out.sorted.bam

rm -f $output_d/$sid.out.sam 
rm -f $output_d/$sid.out.bam

sh $pdir/03-alignment/extract_uniqs.sh   $output_d   $sid.out.sorted.bam


