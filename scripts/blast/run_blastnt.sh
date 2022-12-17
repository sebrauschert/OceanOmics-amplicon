#!/bin/bash

# Blast query to get available taxa

# This is the script to run blastn on the LULU curated DADA2 results using the NCBI nt database
# The query was run seperate for each assay

# Usage: bash scripts/06-blastnt.sh RSV5 16S
voyageID=
assay=
cores=50
#..........................................................................................
usage()
{
          printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay one at a time>\t<string>\n\t\t\t -c <cores, default 50 for blastn>\n\n";
                    exit 1;
}



while getopts v:a:c: flag
do

        case "${flag}" in
                    v) voyageID=${OPTARG};;
                    a) assay=${OPTARG};;
                    c) cores=${OPTARG};;
                    *) usage;;
                esac
done
if [ "${voyageID}" == ""   ]; then usage; fi
if [ "${assay}" == ""   ]; then usage; fi


#activate blast conda environment
eval "$(conda shell.bash hook)"

#"$(conda shell.bash hook)"
conda activate blast-2.12.0

# Now we can use the LULU curated fasta file for the blastn input
echo blasting ${voyageID} ${assay} ASVs

blastn -db /data/tools/databases/ncbi-nt/nt \
       -query 04-LULU/LULU_curated_fasta_${voyageID}_${assay}.fa \
       -num_threads ${cores} \
       -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
       -html > 05-taxa/blast_out/${voyageID}_${assay}_nt.tsv
