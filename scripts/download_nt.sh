#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate blast-2.12.0

# Create database directory in the current working directory
mkdir -p /mnt/scratch/databases/ncbi-nt

cd /mnt/scratch/databases/ncbi-nt

# Download the database. This will take multiple hours
update_blastdb.pl --decompress nt

# Download the taxa database as well
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
# Extract the downloaded files
tar xzvf taxdb.tar.gz
