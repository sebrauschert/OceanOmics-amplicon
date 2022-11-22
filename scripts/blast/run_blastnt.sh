#!/bin/bash

# Blast query to get available taxa

# This is the script to run blastn on the LULU curated DADA2 results using the NCBI nt database
# The query was run seperate for each assay

# Usage: bash scripts/06-blastnt.sh RSV5 16S

#activate amplicon conda environment
eval "$(conda shell.bash hook)"

#"$(conda shell.bash hook)"
conda activate amplicon

# Define the voyage ID and assay in the command line arguments
voyage=$1
assay=$2

#activate blast conda environment
eval "$(conda shell.bash hook)"

#"$(conda shell.bash hook)"
conda activate blast-2.12.0

# Now we can use the LULU curated fasta file for the blastn input
echo blasting ${voyage} ${assay} ASVs

blastn -db /data/tools/databases/ncbi-nt/nt \
       -query 04-LULU/LULU_curated_fasta_${voyage}_${assay}.fa \
       -num_threads 100 \
       -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
       -html > 05-taxa/blast_out/${voyage}_${assay}_nt.tsv
