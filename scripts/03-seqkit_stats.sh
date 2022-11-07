#!/bin/bash

# QC for all data
eval "$(conda shell.bash hook)"
conda activate amplicon

ROOT_DIR=$(pwd)

voyage=$1

# User feedback
echo "Main directory is:"
echo $ROOT_DIR
echo
echo "Statistic files in:"
echo
echo 

# Loop over assays and voyages for rename
for assay in ${@:2} 
    do
    
echo ${ROOT_DIR}/02-demultiplexed/${assay}/

    # Create stats and save to file
    seqkit stats -j 50 -b ${ROOT_DIR}/02-demultiplexed/${assay}/*.fq.gz -a > 01-QC/Sample_statistics_${voyage}_${assay}.txt
 
done
