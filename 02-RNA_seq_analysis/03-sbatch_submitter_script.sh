#!/bin/sh -l
#$ -cwd

. /data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-rlee-v3/richardlee_136891/04-featurecount/config.sh

OUTPUTDIR_FEATURECOUNTS="/data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-rlee-v3/richardlee_136891/04-featurecount"

while read SAMPLEID OUTPUTDIR_HISAT2
do
	mkdir -p $OUTPUTDIR_FEATURECOUNTS/$SAMPLEID

	sbatch --job-name=featureCounts_${SAMPLEID} --output=${SAMPLEID}.slurm.out --error=${SAMPLEID}.std.out  \
           $OUTPUTDIR_FEATURECOUNTS/run_featureCounts.sh \
                   $SAMPLEID        $OUTPUTDIR_FEATURECOUNTS/$SAMPLEID         $OUTPUTDIR_HISAT2

done < RNA_Seq_ids.txt


