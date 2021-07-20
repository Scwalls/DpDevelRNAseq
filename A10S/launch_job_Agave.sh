#!/bin/bash

#SBATCH -n 2
#SBATCH -t 0-12:00
#SBATCH -A scwalls
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

myDir=/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/
sampleFolder=A10S

module load singularity/3.6.3

source /home/scwalls/genome_analysis/GoSTRIPES/bin/xworkStripes -b /scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES -i /scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/gostripes.simg


echo "Launching job"

cd $myDir/$sampleFolder

$rws make -n

$rws make

echo "Job complete"
