#!/bin/bash

#SBATCH -c 1
#SBATCH -t 00-24:00
#SBATCH -o %j.OUT
#SBATCH -e %j.ERROR

module load fastqc/0.11.7

fileDir=/data/LynchLabCME/Daphnia/DaphniaDevel/RNAseq
####### Before running the script, please enter path to desired output directory, below ####
fqDir=/scratch/scwalls/DaphniaDevel/real_RNAseq

echo "Creating links to RNA files."

cd $fqDir

ln -s ${fileDir}/GSF2805-E1_S5_R1_001.fastq.gz E1_S5_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-E1_S5_R2_001.fastq.gz E1_S5_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-E0_S2_R1_001.fastq.gz E0_S2_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-E0_S2_R2_001.fastq.gz E0_S2_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-D1_S4_R1_001.fastq.gz D1_S4_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-D1_S4_R2_001.fastq.gz D1_S4_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-D0_S1_R1_001.fastq.gz D0_S1_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-D0_S1_R2_001.fastq.gz D0_S1_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-C2_S6_R1_001.fastq.gz C2_S6_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-C2_S6_R2_001.fastq.gz C2_S6_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-B1_S3_R1_001.fastq.gz B1_S3_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-B1_S3_R2_001.fastq.gz B1_S3_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-F9_S17_R1_001.fastq.gz F9_S17_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-F9_S17_R2_001.fastq.gz F9_S17_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-F3_S9_R1_001.fastq.gz F3_S9_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-F3_S9_R2_001.fastq.gz F3_S9_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-F10_S20_R1_001.fastq.gz F10_S20_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-F10_S20_R2_001.fastq.gz F10_S20_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-E2_S7_R1_001.fastq.gz E2_S7_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-E2_S7_R2_001.fastq.gz E2_S7_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-D5_S11_R1_001.fastq.gz D5_S11_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-D5_S11_R2_001.fastq.gz D5_S11_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-D4_S10_R1_001.fastq.gz D4_S10_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-D4_S10_R2_001.fastq.gz D4_S10_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-D3_S8_R1_001.fastq.gz D3_S8_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-D3_S8_R2_001.fastq.gz D3_S8_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-C7_S14_R1_001.fastq.gz C7_S14_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-C7_S14_R2_001.fastq.gz C7_S14_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-B9_S16_R1_001.fastq.gz B9_S16_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-B9_S16_R2_001.fastq.gz B9_S16_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-B7_S13_R1_001.fastq.gz B7_S13_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-B7_S13_R2_001.fastq.gz B7_S13_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-B11_S22_R1_001.fastq.gz B11_S22_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-B11_S22_R2_001.fastq.gz B11_S22_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-B10_S19_R1_001.fastq.gz B10_S19_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-B10_S19_R2_001.fastq.gz B10_S19_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-A9_S15_R1_001.fastq.gz A9_S15_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-A9_S15_R2_001.fastq.gz A9_S15_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-A7_S12_R1_001.fastq.gz A7_S12_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-A7_S12_R2_001.fastq.gz A7_S12_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-A1_S23_R1_001.fastq.gz A1_S23_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-A1_S23_R2_001.fastq.gz A1_S23_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-A11_S21_R1_001.fastq.gz A11_S21_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-A11_S21_R2_001.fastq.gz A11_S21_R2_001.fastq.gz
ln -s ${fileDir}/GSF2805-A10_S18_R1_001.fastq.gz A10_S18_R1_001.fastq.gz
ln -s ${fileDir}/GSF2805-A10_S18_R2_001.fastq.gz A10_S18_R2_001.fastq.gz

echo "Starting job"

echo "Running fastqc on the RNA-seq files"

for fq in *.fastq.gz; do
	fastqc $fq
done

echo "Fastqc job is complete"

exit


