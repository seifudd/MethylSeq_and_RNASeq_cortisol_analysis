#!/bin/bash
#$ -cwd

 : <<'END'
inputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter"
outputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/01-results"

cat $inputdir/chr1/RL_chr1_filtered_10.32.cpgs.M.granges.without.SNPs.txt > $outputdir/temp.txt

for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
	cat $outputdir/temp.txt $inputdir/chr$i/RL_chr"$i"_filtered_10.32.cpgs.M.granges.without.SNPs.txt > $outputdir/temp2.txt
	mv -f $outputdir/temp2.txt $outputdir/temp.txt
done

mv -f $outputdir/temp.txt $outputdir/RL_emotional_stress_filtered_10.32.cpgs.M.granges.without.SNPs.txt
rm -f $outputdir/temp1.txt
END

# : <<'END'
inputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter"
outputdir="/home/cluster/nas1/fseifudd/02-methylseq/13-gary-methylseq-human/15-run_bumphunter/01-results"

cat $inputdir/chr1/emotional_stress_high_low_chr1_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt > $outputdir/temp.txt

for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
	awk '{if (NR!=1) {print}}' $inputdir/chr$i/emotional_stress_high_low_chr"$i"_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt > $outputdir/temp1.txt
	cat $outputdir/temp.txt $outputdir/temp1.txt > $outputdir/temp2.txt
	mv -f $outputdir/temp2.txt $outputdir/temp.txt
done

mv -f $outputdir/temp.txt $outputdir/emotional_stress_high_low_bumps_1000.bootstraps_maxgap_300.NOloessbyClustersmoothing.with.covariates.sva.without.snps.zero.constant.added.txt
rm -f $outputdir/temp1.txt
# END



