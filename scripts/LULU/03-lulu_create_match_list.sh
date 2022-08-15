#!/bin/bash

# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate blast-2.12.0

mkdir -p 06-lulu
mkdir -p 06-lulu/database

ln -s $(pwd)/03-dada2/*.fa 06-lulu/database/

cd 06-lulu/database

voyage=$1
declare -a ASSAYS=("16S" "MiFish")

# Loop over assays
for assay in ${ASSAYS[@]} 
  do
  
  #First produce a blastdatabase with the OTUs
  makeblastdb -in ${voyage}_${assay}.fa -parse_seqids -dbtype nucl

  # Then blast the OTUs against the database
  blastn -db ${voyage}_${assay}.fa -outfmt '6 qseqid sseqid pident' -out ../${voyage}_${assay}_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query ${voyage}_${assay}.fa

  # Remove database
  rm -rf *${assay}*.n*
  
done
