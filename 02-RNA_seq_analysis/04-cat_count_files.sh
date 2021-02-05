#!/bin/bash
#$ -cwd


workingdir="/data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-rlee-v3/richardlee_136891/04-featurecount"

for i in `cat RNA_Seq_ids.txt | head -1 | cut -f 1`; do
   echo -e "Gene\t$i" > tmp1
   cat $workingdir/$i/$i.genefeatureCounts.txt | sed '1,2d'  | sort -k1,1 | cut -f 1,7 >> tmp1
done

for i in `cat RNA_Seq_ids.txt | cut -f 1 | sed '1,1d' `; do
   echo -e "$i" > tmp2
   cat $workingdir/$i/$i.genefeatureCounts.txt | sed '1,2d'  | sort -k1,1 | cut -f 7 >> tmp2
   paste tmp1 tmp2 > tmp3
   mv -f tmp3 tmp1
done

for i in `cat RNA_Seq_ids.txt | cut -f 1 `; do
   cat $workingdir/$i/$i.genefeatureCounts.txt | wc -l
done


mv -f $workingdir/tmp1  $workingdir/gene.featurecount.txt
rm -f $workingdir/tmp2
rm -f $workingdir/tmp3

