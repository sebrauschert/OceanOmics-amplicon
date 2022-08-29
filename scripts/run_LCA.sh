#!/bin/bash

# Last Common Ancestor (LCA) analysis to determine most accurate taxa assignments for each ASV
# This is the script to run the CLA script and loop through the voyages, assays and sites

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

	python scripts/07-custom-lca.py \
	--blast_results 04-taxa/blast_out/${voyage}_${assay}_blast_results.tsv \
	--dada2_asv_table 03-dada2/${voyage}/${voyage}_final_table_${assay}.tsv \
	--output 04-taxa/LCA_out/${voyage}_${assay}

 done
done
