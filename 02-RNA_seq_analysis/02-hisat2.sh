

pdir="/data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-rlee-v2/01_batch1_12"

REF="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/grch38/genome"

for sid in `cat $pdir/SIDs` ; do
   sbatch --cpus-per-task=16 --mem=64g  --time=36:00:00  \
      $pdir/03-alignment/sbatch-hisat2.sh   \
      $pdir/01-merged/${sid}.1.fastq.gz  \
      $pdir/01-merged/${sid}.2.fastq.gz  \
      $sid $pdir $REF
done

############# redo for 2 sids

pdir="/data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-rlee-v2/01_batch1_12"

REF="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/grch38/genome"

for sid in 12-0094b 12-0096b ; do
   sbatch --cpus-per-task=16 --mem=64g  --time=36:00:00  \
      $pdir/03-alignment/sbatch-hisat2.sh   \
      $pdir/01-merged/${sid}.1.fastq.gz  \
      $pdir/01-merged/${sid}.2.fastq.gz  \
      $sid $pdir $REF
done


