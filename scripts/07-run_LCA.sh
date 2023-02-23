#!/bin/bash

# Lowest Common Ancestor (LCA) analysis to determine most accurate taxa assignments for each ASV
# This is the script to run the LCA scripts from eDNAFlow and loop through the assays

voyageID=
assay=
database=
#..........................................................................................
usage()
{
          printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay; use flag multiple times for multiple assays>\t<string>\n\t\t\t -d <database; nt, ocom or custom>\t <string>\n\n";
          exit 1;
}
while getopts v:a:d: flag
do

        case "${flag}" in
            v) voyageID=${OPTARG};;
            a) assay+=("$OPTARG");;
            d) database=${OPTARG};;
            *) usage;;
        esac
done
if [ "${voyageID}" == ""  ]; then usage; fi
#if [ "${assay[@]}" == ""  ]; then usage; fi
if [ "${database}" == ""  ]; then usage; fi

#activate blast conda environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

# log the commands
set -x
exec 1>logs/07-run_LCA.log 2>&1

if [ "$database" == "nt" ];
then

for a in ${assay[@]}
  do
        echo  running LCA analysis on ${voyageID} ${a} NCBI nt database

        python scripts/LCA/runAssign_collapsedTaxonomy.py \
        03-dada2/${voyageID}_${a}_lca_input.tsv \
        05-taxa/blast_out/${voyageID}_${a}_nt.tsv \
        100 98 1 \
        05-taxa/LCA_out/${voyageID}_${a}_nt_LCA.tsv
  done

fi


if [ "$database" == "ocom" ];
then

for a in ${assay[@]}
  do
        echo  running LCA analysis on ${voyageID} ${a} OceanOmics database

        python scripts/LCA/runAssign_collapsedTaxonomy.py \
        03-dada2/${voyageID}_${a}_lca_input.tsv \
        05-taxa/blast_out/${voyageID}_${a}_OcOm.tsv \
        100 98 1 \
        05-taxa/LCA_out/${voyageID}_${a}_OcOm_LCA.tsv
  done

fi


if [ "$database" == "custom" ];
then

# load python and taxonkit environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

for a in ${assay[@]}
  do
        echo  running LCA analysis on ${voyageID} ${a} custom database

        python scripts/LCA/runAssign_collapsedTaxonomy.py \
        03-dada2/${voyageID}_${a}_lca_input.tsv \
        05-taxa/blast_out/${voyageID}_${a}_blast_results.tsv \
        100 98 1 \
        05-taxa/LCA_out/${voyageID}_${a}_LCA.tsv
  done

fi
