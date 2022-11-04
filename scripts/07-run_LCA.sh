#!/bin/bash

# Lowest Common Ancestor (LCA) analysis to determine most accurate taxa assignments for each ASV
# This is the script to run the LCA scripts from eDNAFlow and loop through the assays

# Usage: bash scripts/07-run_LCA.sh RSV5 custom 16S MiFish

#activate blast conda environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

# Define voyage ID and database option
voyage=$1
option=$2

if [ "$option" == "nt" ];
then

for assay in ${@:3}
  do
        echo  running LCA analysis on ${voyage} ${assay} NCBI nt database

        python scripts/LCA/runAssign_collapsedTaxonomy.py \
        03-dada2/${voyage}_${assay}_lca_input.tsv \
        05-taxa/blast_out/${voyage}_${assay}_nt.tsv \
        100 98 1 \
        05-taxa/LCA_out/${voyage}_${assay}_nt_LCA.tsv
  done

fi

if [ "$option" == "custom" ];
then

# load python and taxonkit environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

for assay in ${@:3}
  do
        echo  running LCA analysis on ${voyage} ${assay} custom database

        python scripts/LCA/runAssign_collapsedTaxonomy.py \
        03-dada2/${voyage}_${assay}_lca_input.tsv \
        05-taxa/blast_out/${voyage}_${assay}_blast_results.tsv \
        100 98 1 \
        05-taxa/LCA_out/${voyage}_${assay}_LCA.tsv
  done

fi
