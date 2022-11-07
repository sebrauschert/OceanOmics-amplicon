#!/bin/bash

# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate blast-2.12.0

mkdir -p 04-LULU
mkdir -p 04-LULU/database

ln -s $(pwd)/03-dada2/*.fa 04-LULU/database/

cd 04-LULU/database

voyage=$1
assay=$2

#First produce a blastdatabase with the OTUs
makeblastdb -in ${voyage}_${assay}.fa -parse_seqids -dbtype nucl

# Then blast the OTUs against the database
blastn -db ${voyage}_${assay}.fa -outfmt '6 qseqid sseqid pident' -out ../${voyage}_${assay}_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query ${voyage}_${assay}.fa

# Remove database
rm -rf *${assay}*.n*
  
# create lulu asv input table
cd ../
cat ../03-dada2/${voyage}_final_table_${assay}.tsv | tr '\t' ',' > ${voyage}_${assay}_lulu_table.csv
