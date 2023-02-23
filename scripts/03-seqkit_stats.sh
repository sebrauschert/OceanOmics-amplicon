#!/bin/bash

# USAGE:
#   bash 03-seqkit_stats.sh -v <project/voyage ID> \
#                           -a <assay; if multiple assays, then use flag like this: -a 16S -a MiFish> \
#                           -c <number of cores, default 50>

voyageID=
assay=
cores=50

#..........................................................................................
usage()
{
          printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay; use the flag multiple times for multiple assays>\t<string>\n\t\t\t -c <number of cores; default: 50>\n\n";
          exit 1;
}
while getopts v:a:c: flag
do

        case "${flag}" in
            v) voyageID=${OPTARG};;
            a) assay+=("$OPTARG");;
            c) cores=${OPTARG};;
            *) usage;;
        esac
done
if [ "${voyageID}" == ""  ]; then usage; fi
#if [ "${assay}" == ""  ]; then usage; fi

# QC for all data
eval "$(conda shell.bash hook)"
conda activate amplicon

ROOT_DIR=$(pwd)

# log the commands
set -x
echo 'Writing logs to logs/03-seqkit_stats.log'
exec 1>logs/03-seqkit_stats.log 2>&1
# User feedback
echo "Main directory is:"
echo $ROOT_DIR

# Loop over assays and voyages for rename
for a in ${assay[@]}
    do
    
echo ${ROOT_DIR}/01-demultiplexed/${a}/

    # Create stats and save to file
    seqkit stats -j ${cores} -b ${ROOT_DIR}/01-demultiplexed/${a}/*.fq.gz -a > 02-QC/Sample_statistics_${voyageID}_${a}.txt
 
done
