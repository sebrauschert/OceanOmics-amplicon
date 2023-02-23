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

# log the commands
set -x
exec 1>logs/02-rename_demux.log 2>&1

# Rename all demultiplexed files to the sample names

ROOT_DIR=$(pwd)


# User feedback
echo "Main directory is:"
echo $ROOT_DIR
echo

# Loop over assays and voyages for rename
for a in "${assay[@]}"
    do
    # get around small bug where a is empty, leading to nonsense commands
    if [[ -z "${a}" ]];
    then
       continue
    fi
    
    cd ${ROOT_DIR}/01-demultiplexed/${a}
    echo $(pwd)
    
    # Here we reference the rename pattern files in the raw data folder, which contains the information that maps the 
    # index ID pairings with the sample IDs
    mmv < ${ROOT_DIR}/00-raw-data/indices/Sample_name_rename_pattern_${voyageID}_${a}.txt -g
   
    # move the unnamed and unknowns into separate folders 
    mkdir -p ${ROOT_DIR}/01-demultiplexed/${a}/unknown ${ROOT_DIR}/01-demultiplexed/${a}/unnamed
    mv ${ROOT_DIR}/01-demultiplexed/${a}/*unknown*.fq.gz ${ROOT_DIR}/01-demultiplexed/${a}/unknown
    mv ${ROOT_DIR}/01-demultiplexed/${a}/${a}-* ${ROOT_DIR}/01-demultiplexed/${a}/unnamed
     
done
