#!/bin/bash

#SBATCH -c 1                        
#SBATCH -t 01-6:00
#SBATCH -o %j.OUT
#SBATCH -e %j.ERROR

module load star/2.4.2a

fileDir=/data/LynchLabCME/Daphnia/DaphniaDevel/RNAseq
####### Before running the script, please enter path to desired output directory, below ####
WD=/scratch/scwalls/DaphniaDevel/STAR_RNA
fqDir=fastqs
genomedir=/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/DpGENOME
genomeFasta=PA42.4.1.fasta
configfile=/scratch/scwalls/T502_RNAseq/scripts/STARalign.conf
nThreads=16


function readconfigfile {
# Read the specified ($1) STARalign configuration file:
if [ ! -e "$1" ] ; then
  echo ""
  echo "Fatal error: STARalign config file $1 does not exist. Please check."
  exit 1
fi

STARgenomeGenerateOptions=`grep '^STARgenomeGenerateOptions=' "$1" | awk -F"=" '{print $2}'`
STARalignReadsOptions=`grep '^STARalignReadsOptions=' "$1" | awk -F"=" '{print $2}'`
bamFilterOptions=`grep '^bamFilterOptions=' "$1" | awk -F"=" '{print $2}'`
}

echo "Starting job"

echo "Making symbolic links to fastq files."

cd $WD

if [ ! -d "$fqDir" ] ; then
    mkdir $fqDir
  fi

cd $fqDir

ln -s ${fileDir}/GSF2805-E1_S5_R1_001.fastq.gz E1_S5_R1.fastq.gz
ln -s ${fileDir}/GSF2805-E1_S5_R2_001.fastq.gz E1_S5_R2.fastq.gz
ln -s ${fileDir}/GSF2805-E0_S2_R1_001.fastq.gz E0_S2_R1.fastq.gz
ln -s ${fileDir}/GSF2805-E0_S2_R2_001.fastq.gz E0_S2_R2.fastq.gz
ln -s ${fileDir}/GSF2805-D1_S4_R1_001.fastq.gz D1_S4_R1.fastq.gz
ln -s ${fileDir}/GSF2805-D1_S4_R2_001.fastq.gz D1_S4_R2.fastq.gz
ln -s ${fileDir}/GSF2805-D0_S1_R1_001.fastq.gz D0_S1_R1.fastq.gz
ln -s ${fileDir}/GSF2805-D0_S1_R2_001.fastq.gz D0_S1_R2.fastq.gz
ln -s ${fileDir}/GSF2805-C2_S6_R1_001.fastq.gz C2_S6_R1.fastq.gz
ln -s ${fileDir}/GSF2805-C2_S6_R2_001.fastq.gz C2_S6_R2.fastq.gz
ln -s ${fileDir}/GSF2805-B1_S3_R1_001.fastq.gz B1_S3_R1.fastq.gz
ln -s ${fileDir}/GSF2805-B1_S3_R2_001.fastq.gz B1_S3_R2.fastq.gz
ln -s ${fileDir}/GSF2805-F9_S17_R1_001.fastq.gz F9_S17_R1.fastq.gz
ln -s ${fileDir}/GSF2805-F9_S17_R2_001.fastq.gz F9_S17_R2.fastq.gz
ln -s ${fileDir}/GSF2805-F3_S9_R1_001.fastq.gz F3_S9_R1.fastq.gz
ln -s ${fileDir}/GSF2805-F3_S9_R2_001.fastq.gz F3_S9_R2.fastq.gz
ln -s ${fileDir}/GSF2805-F10_S20_R1_001.fastq.gz F10_S20_R1.fastq.gz
ln -s ${fileDir}/GSF2805-F10_S20_R2_001.fastq.gz F10_S20_R2.fastq.gz
ln -s ${fileDir}/GSF2805-E2_S7_R1_001.fastq.gz E2_S7_R1.fastq.gz
ln -s ${fileDir}/GSF2805-E2_S7_R2_001.fastq.gz E2_S7_R2.fastq.gz
ln -s ${fileDir}/GSF2805-D5_S11_R1_001.fastq.gz D5_S11_R1.fastq.gz
ln -s ${fileDir}/GSF2805-D5_S11_R2_001.fastq.gz D5_S11_R2.fastq.gz
ln -s ${fileDir}/GSF2805-D4_S10_R1_001.fastq.gz D4_S10_R1.fastq.gz
ln -s ${fileDir}/GSF2805-D4_S10_R2_001.fastq.gz D4_S10_R2.fastq.gz
ln -s ${fileDir}/GSF2805-D3_S8_R1_001.fastq.gz D3_S8_R1.fastq.gz
ln -s ${fileDir}/GSF2805-D3_S8_R2_001.fastq.gz D3_S8_R2.fastq.gz
ln -s ${fileDir}/GSF2805-C7_S14_R1_001.fastq.gz C7_S14_R1.fastq.gz
ln -s ${fileDir}/GSF2805-C7_S14_R2_001.fastq.gz C7_S14_R2.fastq.gz
ln -s ${fileDir}/GSF2805-B9_S16_R1_001.fastq.gz B9_S16_R1.fastq.gz
ln -s ${fileDir}/GSF2805-B9_S16_R2_001.fastq.gz B9_S16_R2.fastq.gz
ln -s ${fileDir}/GSF2805-B7_S13_R1_001.fastq.gz B7_S13_R1.fastq.gz
ln -s ${fileDir}/GSF2805-B7_S13_R2_001.fastq.gz B7_S13_R2.fastq.gz
ln -s ${fileDir}/GSF2805-B11_S22_R1_001.fastq.gz B11_S22_R1.fastq.gz
ln -s ${fileDir}/GSF2805-B11_S22_R2_001.fastq.gz B11_S22_R2.fastq.gz
ln -s ${fileDir}/GSF2805-B10_S19_R1_001.fastq.gz B10_S19_R1.fastq.gz
ln -s ${fileDir}/GSF2805-B10_S19_R2_001.fastq.gz B10_S19_R2.fastq.gz
ln -s ${fileDir}/GSF2805-A9_S15_R1_001.fastq.gz A9_S15_R1.fastq.gz
ln -s ${fileDir}/GSF2805-A9_S15_R2_001.fastq.gz A9_S15_R2.fastq.gz
ln -s ${fileDir}/GSF2805-A7_S12_R1_001.fastq.gz A7_S12_R1.fastq.gz
ln -s ${fileDir}/GSF2805-A7_S12_R2_001.fastq.gz A7_S12_R2.fastq.gz
ln -s ${fileDir}/GSF2805-A1_S23_R1_001.fastq.gz A1_S23_R1.fastq.gz
ln -s ${fileDir}/GSF2805-A1_S23_R2_001.fastq.gz A1_S23_R2.fastq.gz
ln -s ${fileDir}/GSF2805-A11_S21_R1_001.fastq.gz A11_S21_R1.fastq.gz
ln -s ${fileDir}/GSF2805-A11_S21_R2_001.fastq.gz A11_S21_R2.fastq.gz
ln -s ${fileDir}/GSF2805-A10_S18_R1_001.fastq.gz A10_S18_R1.fastq.gz
ln -s ${fileDir}/GSF2805-A10_S18_R2_001.fastq.gz A10_S18_R2.fastq.gz

cd $WD

echo "Retrieving Daphnia  genome assembly and annotation files"

#source 0README

echo "Performing alignment on the RNA-seq files"

readconfigfile $configfile

cd $genomedir

  if [ ! -e ${genomedir}/SAindex ] ; then
    echo "STAR --runMode genomeGenerate --runThreadN $numproc  ${STARgenomeGenerateOptions}  --genomeDir $genomedir --genomeFastaFiles $genomedir/*.fa"
    STAR --runMode genomeGenerate --runThreadN $nThreads  ${STARgenomeGenerateOptions}  --genomeDir $genomedir --genomeFastaFiles $genomedir/${genomeFasta}
  else
    echo "Using existing STAR suffix array for genome file $genomedir/*.fa"
  fi

  cd ${WD}/${fqDir}
  for file1 in *_R1.fastq.gz; do
    file2=$(basename $file1 _R1.fastq.gz)_R2.fastq.gz
    echo "STAR --runMode alignReads --runThreadN $numproc  ${STARalignReadsOptions}  --outSAMtype BAM SortedByCoordinate --outSAMorder Paired  --outFileNamePrefix $(basename $file1 _001.fastq.gz).STAR.  --genomeDir $genomedir  --readFilesIn ${file1} ${file2}"
    STAR --runMode alignReads --runThreadN $nThreads  ${STARalignReadsOptions}  --outSAMtype BAM SortedByCoordinate --outSAMorder Paired  --outFileNamePrefix $(basename $file1 .fastq).STAR. --genomeDir $genomedir  --readFilesIn ${file1} ${file2}
  done

  cd ${WD}
  if [ ! -d "alignments" ] ; then
    mkdir alignments
  fi
  cd alignments
  for file1 in ${WD}/$fqDir/*.STAR.*.bam ; do
    echo $file1
    file2=$(basename $file1)
    echo $file2
    ln -s $file1 ./${file2/.STAR.Aligned.sortedByCoord.out}
  done
  cd ..

  echo ""
  echo " Done with step 3 (read mapping)."
  echo ""
  echo "================================================================================"
fi


exit


