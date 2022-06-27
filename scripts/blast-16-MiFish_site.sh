#!/bin/bash

# Blast query to get available taxa using the specific 16S and MiFish databases
# This script loops through all voyages, assays and sites for each sequencing run, calling on the blast script: 06-blast-16S-MiFish.py

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

	declare -a SITES=$(ls 02-demultiplexed/${voyage}/${assay} | grep -v "gz" | grep -v "unknown" | grep -v "unnamed" | grep -v "README.md")
	echo  02-demultiplexed/${voyage}/${assay}
	echo $SITES

	for site in ${SITES[@]}
	 do

python scripts/06-blast-16S-MiFish.py \
	 --dada2_file 03-dada2/${voyage}/${voyage}_${assay}_${site}.fa \
	 --out_path 04-taxa/blast_out/${voyage}_${site}_ \
	 --database ${assay}

  done
 done
done
