#!/bin/bash

# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate blast-2.12.0

mkdir -p 06-lulu
mkdir -p 06-lulu/database

ln -s $(pwd)/03-dada2/*.fa 06-lulu/database/

cd 06-lulu/database

#First produce a blastdatabase with the OTUs
makeblastdb -in RS_16S_merged.fa -parse_seqids -dbtype nucl

# Then blast the OTUs against the database
blastn -db RS_16S_merged.fa -outfmt '6 qseqid sseqid pident' -out ../RS_16S_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query RS_16S_merged.fa

# Remove database
rm -rf *16S*.n*

#First produce a blastdatabase with the OTUs
makeblastdb -in RS_MiFish_merged.fa -parse_seqids -dbtype nucl

# Then blast the OTUs against the database
blastn -db RS_MiFish_merged.fa -outfmt '6 qseqid sseqid pident' -out ../RS_MiFish_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query RS_MiFish_merged.fa


# Remove database
rm -rf *MiFish*.n*