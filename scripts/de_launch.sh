#!/bin/bash

#SBATCH -c 1
#SBATCH -t 01-12:00
#SBATCH -o %j.OUT
#SBATCH -e %j.ERROR

module load r/4.0.2

WD=/scratch/scwalls/T502_RNAseq/scripts
#WD=/N/dc2/scratch/rtraborn/T502_RNAseq/scripts

cd $WD

echo "Launching DE expression job"

R CMD BATCH de_analysis_Daphnia_AB.R

echo "Done."

exit
