#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-4:00
#SBATCH -A scwalls
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load r/4.0.2

myDir=/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/tsr_Es

echo "Launching job"

cd $myDir

./xrunSwf > err

echo "Job complete"
