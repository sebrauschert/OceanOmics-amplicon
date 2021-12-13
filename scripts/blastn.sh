#!/bin/bash

# Blast query to get available taxa

# This is the script to run blast on the DADA2 results
# It was run on my nimbus instance and the nt database was downloaded into
# /data/ubuntu/databases.
# The query was run seperate by fish and elasmo


blastn -db databases/nt \
       -query <query.fa> \
       -evalue 1e-3 \
       -num_threads <number of threads to be used> \
       -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
       -html > tabular_output.tsv
