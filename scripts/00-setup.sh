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

# from https://stackoverflow.com/a/246128
# where does the setup script live?
# need to know this for more robust copying of all scripts later in the script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Create a new directory based on input
#...............................................................................................
mkdir ${projectID}_amplicon_analysis

# Enter the folder and make it a datalad container
cd ${projectID}_amplicon_analysis
pwd
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
      01-demultiplexed \
      02-QC \
      03-dada2/QC_plots \
      03-dada2/tmpfiles \
      03-dada2/errorModel \
      04-LULU \
      05-taxa/blast_out \
      05-taxa/LCA_out \
      06-report \
      databases \
      scripts \
      logs

# Place a README.md in every folder
#...............................................................................................
touch README.md

echo "# Step:" >> README.md
echo "# Analyst:" >> README.md
echo "# Data locations:" >> README.md
echo "# Script used:" >> README.md
echo "# Software version:" >> README.md
echo "# Problems encountered:" >> README.md


parallel cp README.md ::: 01-demultiplexed \
      02-QC \
      03-dada2 \
      04-LULU \
      05-taxa \
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

cp -r ${SCRIPT_DIR}/* scripts

# Finished
#...............................................................................................
echo "Finished setting up analysis directory!"
echo
echo "Directory name is `pwd`"
echo
tree -d
#...............................................................................................
