#!/bin/bash

# Blast query to get available taxa

# This is the script to run blast on the DADA2 results
# It was run on my nimbus instance and the nt database was downloaded into
# /data/ubuntu/databases.
# The query was run seperate for 16S and MiFish assays

#activate blast conda environment
source /data/software/miniconda3/etc/profile.d/conda.sh

#"$(conda shell.bash hook)"
conda activate blast-2.12.0

blastn -db /DATA/databases/nt \
       -query $1 \
       -evalue 1e-3 \
       -num_threads $2 \
       -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
       -html > 06-taxa/tabular_output_$3.tsv
