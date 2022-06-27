#!/bin/bash

# Lowest Common Ancestor (LCA) analysis to determine most accurate taxa assignments for each ASV
# This script loops through all voyages, assays and sites for this sequencing run, calling on the LCA script: 07-custom-lca.py

#activate blast conda environment
source /home/jpearce/data/tools/miniconda3_jpearce/etc/profile.d/conda.sh

#"$(conda shell.bash hook)"
conda activate pytaxonkit

# name voyage and assay
declare -a VOYAGES=$(ls 02-demultiplexed | grep -v "README.md")
declare -a ASSAYS=("16S" "MiFish")

for voyage in ${VOYAGES[@]}
  do

  for assay in ${ASSAYS[@]}
    do
declare -a SITES=$(ls 02-demultiplexed/${voyage}/${assay} | grep -v "gz" | grep -v "unknown" | grep -v "unnamed" | grep -v README)

echo $SITES
echo  02-demultiplexed/${voyage}/${assay}

        for site in ${SITES[@]}
         do

	python scripts/07-custom-lca.py \
	--blast_results 04-taxa/blast_out/${site}_${assay}_blast_results.tsv \
	--dada2_asv_table 03-dada2/${voyage}/${voyage}_asv_final_table_${assay}_${site}.tsv \
	--output 04-taxa/LCA_out/${voyage}_${assay}_${site}

  done
 done
done
