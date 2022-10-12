#!/bin/bash

# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

voyage=$1
assay=$2

# Usage: bash scripts/LULU/01-get_lineage.sh PCV3 MiFish

# Provide as input the filename with the taxid

cat 04-taxa/${voyage}_${assay}_decontam_table.csv \
| tr ',' '\t' \
| cut -f 2 \
|  awk 'NR>1' \
|  taxonkit lineage \
| taxonkit reformat \
| csvtk -H -t cut -f 1,3 \
|  csvtk -H -t sep -f 2 -s ';' -R \
|  csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species > 04-taxa/taxonkit_lineage_${voyage}_${assay}.tsv