#!/bin/bash

# Blast query to get available taxa

# This is the script to run blastn on the LULU curated DADA2 results using the NCBI nt database
# The query was run seperate for each assay

voyageID=
#assay=
database=
cores=50
#..........................................................................................
usage()
{
          printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay; use flag multiple times for multiple assays>\t<string>\n\t\t\t -d <database; nt, ocom or custom>\t <string>\n\t\t\t -c <cores, default 50 for blastn>\n\n";
          exit 1;
}
while getopts v:a:d:c: flag
do

        case "${flag}" in
            v) voyageID=${OPTARG};;
            a) assay+=("$OPTARG");;
            d) database=${OPTARG};;
            c) cores=${OPTARG};;
            *) usage;;
        esac
done
if [ "${voyageID}" == ""  ]; then usage; fi
#if [ "${assay}" == ""  ]; then usage; fi
if [ "${database}" == ""  ]; then usage; fi

# log the commands
set -x
exec 1>logs/06-run_blast.log 2>&1

# Define voyage ID and database option

if [ "${database}" == "nt" ];
then

for a in ${assay[@]}
  do
  bash scripts/blast/run_blastnt.sh -v ${voyageID} -a ${a} -c ${cores}
  done

fi


if [ "${database}" == "ocom" ];
then

for a in ${assay[@]}
  do
  bash scripts/blast/run_blastOcOm.sh -v ${voyageID} -a ${a} -c ${cores}
  done

fi


if [ "${database}" == "custom" ];
then

# load python and taxonkit environment
eval "$(conda shell.bash hook)"
conda activate pytaxonkit

for a in ${assay[@]}
  do
  python scripts/blast/blast-16S-MiFish.py \
         --dada2_file 04-LULU/LULU_curated_fasta_${voyageID}_${a}.fa \
         --out_path 05-taxa/blast_out/${voyageID}_ \
         --database ${a}
  done

fi
