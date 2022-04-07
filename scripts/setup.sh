#!/bin/bash

# Source conda so we can use it with environments
#...............................................................................................
eval "$(conda shell.bash hook)"

# 1.) Create the directory structure; for this we require the datalad software to be installed 
# https://www.datalad.org/
#...............................................................................................

conda activate datalad

# Create a new directory based on input
#...............................................................................................
mkdir $(basename $1)_Amplicon_$(whoami)

# Enter the folder and make it a datalad container
cd $(basename $1)_Amplicon_$(whoami)

echo 'Preparing data repository...'
echo ''

datalad create .

# And finally, add the necessary lines to the .gitattributes that we need
#...............................................................................................
echo -en '\n***.org annex.largefiles=nothing' >> .gitattributes
echo -en '\n***.sh annex.largefiles=nothing' >> .gitattributes
echo -en '\n***.txt annex.largefiles=nothing' >> .gitattributes
echo -en '\n***.py annex.largefiles=nothing' >> .gitattributes
echo -en '\n***.ipynb annex.largefiles=nothing' >> .gitattributes
echo -en '\n***.R annex.largefiles=nothing' >> .gitattributes
echo -en '\n***.Rmd annex.largefiles=nothing' >> .gitattributes
echo -en '\n***.md annex.largefiles=nothing' >> .gitattributes
echo -en '\n** annex.largefiles(largerthan=200kb)' >> .gitattributes

# We need to save the changes now
#...............................................................................................   
datalad save -m "Appended .gitattributes to not track text files"
    
# Finished
echo ''
echo ''
echo 'data repository created...'
echo 'path is:' $PWD
    

# Set up the directory structure
#...............................................................................................
mkdir -p 00-raw-data/indices \
      01-QC \
      02-demultiplexed/sample_names \
      03-dada2/QC_plots \
      04-taxa \
      05-report \
      scripts

# Place a README.md in every folder
#...............................................................................................
touch README.md

echo "# step:" >> README.md
echo "# analyst" >> README.md
echo "# data locations" >> README.md
echo "# script used" >> README.md
echo "# software version" >> README.md
echo "# problems encountered" >> README.md

conda deactivate

parallel cp README.md ::: 01-QC \
      02-demultiplexed \
      03-dada2 \
      04-taxa \
      05-report

# Remove the readme file from the main folder structure
#...............................................................................................
rm README.md

# Create a general README for this project
#...............................................................................................
touch README.md
echo "# project" >> README.md
echo "# analyst" >> README.md
echo "# overview" >> README.md

cp -r ../OceanOmics-amplicon/scripts .

conda activate datalad

datalad save -m "Analysis repo setup"

# Finished
#...............................................................................................
echo "Finished setting up analysis directory!\n"
echo "Directory name is $1"
echo
tree -d
#...............................................................................................