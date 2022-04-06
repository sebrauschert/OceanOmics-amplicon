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
    mv 16S-*16S-*.fq.gz unnamed
    
  
  done
done

# Error message RS21 16S
#16S-55F-16S-60R.R1.fq.gz , 16S-55F-16S-61R.R1.fq.gz -> RS-CL-BL-3_16S.1.fq.gz : collision.
#16S-55F-16S-60R.R2.fq.gz , 16S-55F-16S-61R.R2.fq.gz -> RS-CL-BL-3_16S.2.fq.gz : collision.
#16S-56F-16S-59R.R1.fq.gz , 16S-56F-16S-64R.R1.fq.gz -> RS-CL-S2-4_16S.1.fq.gz : collision.
#16S-56F-16S-59R.R2.fq.gz , 16S-56F-16S-64R.R2.fq.gz -> RS-CL-S2-4_16S.2.fq.gz : collision.

# Error message RS21 MiFish
# MiFish-F15-MiFish-R4.R1.fq.gz , MiFish-F15-MiFish-R5.R1.fq.gz -> RS-CL-BL-3_MiFish.1.fq.gz : collision.
# MiFish-F15-MiFish-R4.R2.fq.gz , MiFish-F15-MiFish-R5.R2.fq.gz -> RS-CL-BL-3_MiFish.2.fq.gz : collision.
# MiFish-F16-MiFish-R3.R1.fq.gz , MiFish-F16-MiFish-R8.R1.fq.gz -> RS-CL-S2-4_MiFish.1.fq.gz : collision.
# MiFish-F16-MiFish-R3.R2.fq.gz , MiFish-F16-MiFish-R8.R2.fq.gz -> RS-CL-S2-4_MiFish.2.fq.gz : collision.

