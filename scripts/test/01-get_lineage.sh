#!/bin/bash


# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

# Provide as input the filename with the taxid

cat $1 \
| tr ',' '\t' \
| cut -f 2 \
| head \
|  awk 'NR>1' \
|  taxonkit lineage \
| taxonkit reformat \
| csvtk -H -t cut -f 1,3 \
|  csvtk -H -t sep -f 2 -s ';' -R \
|  csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species > 04-taxa/taxonkit_lineage_$(basename $1 .csv).tsv