#!/bin/bash

# QC for all data
eval "$(conda shell.bash hook)"
conda activate amplicon

ROOT_DIR=$(pwd)

declare -a VOYAGES=$(ls 02-demultiplexed | grep -v "README.md" | grep -v "sample_names")
declare -a ASSAYS=("16S" "MiFish")

# User feedback
echo "Main directory is:"
echo $ROOT_DIR
echo
echo "Statistic files in:"
echo


# Loop over assays and voyages for rename
for voyage in ${VOYAGES[@]}
  do

  for assay in ${ASSAYS[@]} 
    do
    
    # Create stats and save to file
    seqkit stats -j 50 -b ${ROOT_DIR}/02-demultiplexed/${voyage}/${assay}/*.fq.gz -a > 01-QC/Sample_statistics_${voyage}_${assay}.txt
 
  done
done