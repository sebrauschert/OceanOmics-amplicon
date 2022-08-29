#!/bin/bash

# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

# Usage: bash scripts/08-LULU.sh

# name voyage and assay
declare -a VOYAGES=$(ls 02-demultiplexed | grep -v "README.md")
declare -a ASSAYS=("16S" "MiFish")

for voyage in ${VOYAGES[@]}
  do

  for assay in ${ASSAYS[@]}
    do
        echo  ${voyage}/${assay}
        
        bash scripts/LULU/01-get_lineage.sh ${voyage} ${assay}
        
        Rscript scripts/LULU/02-merge_lineage_with_LCA.R -v ${voyage} -a ${assay}
        
        bash scripts/LULU/03-lulu_create_match_list.sh ${voyage}
        
        Rscript scripts/LULU/04-LULU.R -v ${voyage} -a ${assay}
        
        Rscript scripts/LULU/05-create_phyloseq_object.R -v ${voyage} -a ${assay}

 done
done
