#!/bin/sh -l
#$ -cwd

#SBATCH --job-name="wake_runbumphunter"
#SBATCH --partition=largemem
#SBATCH --time=72:00:00
#SBATCH --mem=1000g
#SBATCH --cpus-per-task=1

scriptdir="/data/NHLBI_BCB/Fayaz/19-rales12-bumphunter-results-dichot"

cd $scriptdir
module load R
R CMD BATCH ./bumphunter.wake.tert.marcc.R

