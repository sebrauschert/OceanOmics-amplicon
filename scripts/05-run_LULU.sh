#!/bin/bash

# load python and taxonkit environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

# Usage: bash scripts/05-run_LULU.sh RSV5 16S MiFish

# name voyage and assay
voyage=$1

for assay in ${@:2}
  do
      echo  "Running LULU on ${voyage} ${assay}"
        
        bash scripts/LULU/01-lulu_create_match_list.sh ${voyage} ${assay}
        
        Rscript scripts/LULU/02-LULU.R -v ${voyage} -a ${assay}
        
        # Next we need to curate the fasta files from DADA2 to only include the ASVs output by LULU
        echo curating ${voyage} ${assay} fasta file

        cat 04-LULU/LULU_curated_counts_${voyage}_${assay}.csv | \
        cut -d "," -f1 | \
        sed 1,1d | \
        seqkit grep -f - 03-dada2/${voyage}_${assay}.fa -o 04-LULU/LULU_curated_fasta_${voyage}_${assay}.fa

done
