#!/bin/bash

# Blast query to get available taxa

# This is the script to run blastn on the LULU curated DADA2 results using the NCBI nt database
# The query was run seperate for each assay

# Usage: bash scripts/06-blastnt.sh RSV5 custom 16S MiFish

# Define voyage ID and database option
voyage=$1
option=$2

if [ "$option" == "nt" ];
then

for assay in ${@:3}
  do
  bash scripts/run_blastnt.sh ${voyage} ${assay}
  done

fi

if [ "$option" == "custom" ];
then

# load python and taxonkit environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

for assay in ${@:3}
  do
  python scripts/blast-16S-MiFish.py \
         --dada2_file 04-LULU/LULU_curated_fasta_${voyage}_${assay}.fa \
         --out_path 05-taxa/blast_out/${voyage}_ \
         --database ${assay}
  done

fi
