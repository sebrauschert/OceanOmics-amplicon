#!/bin/bash

# Rename all demultiplexed files to the sample names

ROOT_DIR=$(pwd)

declare -a VOYAGES=$(ls 02-demultiplexed | grep -v "README.md" | grep -v "sample_names")
declare -a ASSAYS=("16S" "MiFish")

# User feedback
echo "Main directory is:"
echo $ROOT_DIR
echo
echo "Rename files in:"
echo

# Loop over assays and voyages for rename
for voyage in ${VOYAGES[@]}
  do

  for assay in ${ASSAYS[@]} 
    do
    
    cd ${ROOT_DIR}/02-demultiplexed/${voyage}/${assay}
    echo $(pwd)
    
    # Here we reference the rename pattern files in the raw data folder, which contains the information that maps the 
    # index ID pairings with the sample IDs
    mmv < ${ROOT_DIR}/02-demultiplexed/sample_names/Sample_name_rename_pattern_${voyage}_${assay}.txt
   
    # move the unnamed and unknowns into separate folders 
    mkdir unknown unnamed
    mv *unknown*.fq.gz unknown
    mv ${assay}-*${assay}-*.fq.gz unnamed
    
  
  done
done