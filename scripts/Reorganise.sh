#!/bin/bash

# Put control samples into a seperate Controls folder

#..........................................................................................
usage()
{
          printf "Usage: $0 -a <assay; use the flag multiple times for multiple assays>\t<string>\n\n";
          exit 1;
}
while getopts a: flag
do

        case "${flag}" in
            a) assay+=("$OPTARG");;
            *) usage;;
        esac
done

for a in "${assay[@]}"
do
    if [[ -z "${a}" ]];
    then
       continue
    fi

    input_directory="$(pwd)/01-demultiplexed/${a}"

    mkdir ${input_directory}/Controls
    mv ${input_directory}/*EB* ${input_directory}/Controls
    mv ${input_directory}/*WC* ${input_directory}/Controls
    mv ${input_directory}/*FC* ${input_directory}/Controls
done
