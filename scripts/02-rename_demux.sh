#!/bin/bash



voyageID=
assay=
#..........................................................................................
usage()
{
          printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay; use flag multiple times for multiple assays>\t<string>\n\n";
          exit 1;
}
while getopts v:a: flag
do

        case "${flag}" in
            v) voyageID=${OPTARG};;
            a) assay+=("$OPTARG");;
            *) usage;;
        esac
done
if [ "${voyageID}" == ""  ]; then usage; fi
#if [ "${assay}" == ""  ]; then usage; fi

# Rename all demultiplexed files to the sample names

ROOT_DIR=$(pwd)

echo ${voyageID}

# User feedback
echo "Main directory is:"
echo $ROOT_DIR
echo

# Loop over assays and voyages for rename
for a in "${assay[@]}"
    do
    
    cd ${ROOT_DIR}/01-demultiplexed/${a}
    echo $(pwd)
    
    # Here we reference the rename pattern files in the raw data folder, which contains the information that maps the 
    # index ID pairings with the sample IDs
    mmv < ${ROOT_DIR}/00-raw-data/indices/Sample_name_rename_pattern_${voyageID}_${a}.txt
   
    # move the unnamed and unknowns into separate folders 
    mkdir -p unknown unnamed
    mv *unknown*.fq.gz unknown
    mv ${a}-* unnamed
     
done
