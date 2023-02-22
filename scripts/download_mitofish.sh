#!/usr/bin/env bash

# Create database directory in the current working directory
mkdir -p /mnt/scratch/databases/MiFish

cd /mnt/scratch/databases/MiFishDB

# Download and fix
wget http://mitofish.aori.u-tokyo.ac.jp/files/mitogenomes.zip
unzip mitogenomes.zip
rm *genes.fa
cat *.fa > MiFishDB.fasta
