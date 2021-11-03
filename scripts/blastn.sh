#!/bin/bash

# Blast query to get available taxa

blastn -db databases/nt \
  -query voyage1elasmo.fa \
  -evalue 1e-3 \
  -num_threads 16 \
  -outfmt "6 qseqid sseqid sacc stitle sscinames staxids ssciname scomname sskingdoms sblastnames pident slen length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -html > voyage1_elasmo_results_2_tabular