#!/bin/bash
set -e
set -o pipefail # see http://redsymbol.net/articles/unofficial-bash-strict-mode/
# set -u does not work well with SLURM, conda, or HPC module systems IME

# USAGE:
# bash 00-setup.sh -v <project/voyage ID>

projectID=

usage()
{
          printf "Usage: $0 -p <projectID>\t<string>\n\n";
          exit 1;
}
while getopts p: flag
do

        case "${flag}" in
            p) projectID=${OPTARG};;
            *) usage;;
        esac
done
if [ "${projectID}" == ""  ]; then usage; fi


# Create a new directory based on input
#...............................................................................................
mkdir mnt/scratch/${projectID}_amplicon_analysis

# Enter the folder and make it a datalad container
cd /mnt/scratch/${projectID}_amplicon_analysis

echo 'Preparing data repository...'
echo ''

# Finished
echo ''
echo ''
echo 'data repository created...'
echo 'path is:' $(pwd)


# Set up the directory structure
#...............................................................................................
mkdir -p 00-raw-data/indices \
      01-demultiplexed/$2 \
      01-demultiplexed/$3 \
      02-QC \
      03-dada2/QC_plots \
      03-dada2/tmpfiles \
      03-dada2/errorModel \
      04-taxa/blast_out \
      04-taxa/LCA_out \
      05-LULU \
      06-report \
      scripts

# Place a README.md in every folder
#...............................................................................................
touch README.md

echo "# Step:" >> README.md
echo "# Analyst:" >> README.md
echo "# Data locations:" >> README.md
echo "# Script used:" >> README.md
echo "# Software version:" >> README.md
echo "# Problems encountered:" >> README.md


parallel cp README.md ::: 02-QC \
      01-demultiplexed \
      03-dada2 \
      04-taxa \
      05-LULU \
      06-report

# Remove the readme file from the main folder structure
#...............................................................................................
rm README.md

# Create a general README for this project
#...............................................................................................
touch README.md
echo "# Project:" >> README.md
echo "# Analyst:" >> README.md
echo "# Overview:" >> README.md


# Removing the copying - we assume that we're inside the OceanOmics-amplicon folder anyways
# get current folder
cp -r ../OceanOmics-amplicon/scripts/* scripts

# Finished
#...............................................................................................
echo "Finished setting up analysis directory!\n"
echo "Directory name is $1"
echo
tree -d
#...............................................................................................
