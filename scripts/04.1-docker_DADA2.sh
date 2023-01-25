#!/bin/bash

# USAGE
# bash DADA2.sh -v <project/voyage ID> \
#               -a <assay; use the flag multiple timed for multiple assays, e.g. -a 16S -a MiFish> \
#               -o <option; pooled> \
#               -c <number of cores, default 50>

voyageID=
#assay=
option=
cores=50

#..........................................................................................
usage()
{
          printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay; use the flag multiple times for multiple assays>\t<string>\n\t\t\t -c <number of cores; default: 50>\n\n";
          exit 1;
}
while getopts v:a:o:c: flag
do

        case "${flag}" in
            v) voyageID=${OPTARG};;
            a) assay+=("$OPTARG");;
            o) option=${OPTARG};;
            c) cores=${OPTARG};;
            *) usage;;
        esac
done
if [ "${voyageID}" == ""  ]; then usage; fi
#if [ "${assay[@]}" == ""  ]; then usage; fi

Rscript /opt/amplicon_pipeline/04-DADA2.R -v $voyageID -a $assay -o $option -c $cores
