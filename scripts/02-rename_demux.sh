#!/bin/bash

# Rename all demultiplexed files to the sample names

ROOT_DIR=$(pwd)

voyage=$1
echo "$1"

# User feedback
echo "Main directory is:"
echo $ROOT_DIR
echo
echo "Rename files in:"
echo 

# Loop over assays and voyages for rename
for assay in ${@:2} 
    do
    
    cd ${ROOT_DIR}/01-demultiplexed/${assay}
    echo $(pwd)
    
    # Here we reference the rename pattern files in the raw data folder, which contains the information that maps the 
    # index ID pairings with the sample IDs
    mmv < ${ROOT_DIR}/00-raw-data/indices/Sample_name_rename_pattern_${voyage}_${assay}.txt
   
    # move the unnamed and unknowns into separate folders 
    mkdir -p unknown unnamed
    mv *unknown*.fq.gz unknown
    mv ${assay}-* unnamed
     
done
