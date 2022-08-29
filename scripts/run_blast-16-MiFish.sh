#!/bin/bash

# Blast query to get available taxa

# This is the script to run blast on the DADA2 results
# It was run on my nimbus instance and the nt database was downloaded into
# /data/ubuntu/databases.
# The query was run seperate for 16S and MiFish assays

#activate blast conda environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

# name voyage and assay
declare -a VOYAGES=("PCV3")
declare -a ASSAYS=("MiFish")

for voyage in ${VOYAGES[@]}
  do

  for assay in ${ASSAYS[@]}
    do

	echo  ${voyage}/${assay}

python scripts/06-blast-16S-MiFish.py \
	 --dada2_file 03-dada2/${voyage}/${voyage}_${assay}.fa \
	 --out_path 04-taxa/blast_out/${voyage}_ \
	 --database ${assay}

  done
done
