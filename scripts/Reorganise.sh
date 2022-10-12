#!/bin/bash

# Organising samples into sites fro DADA2 input
# The query was run seperate for 16S and MiFish assays

# name voyage and assay
voyage=$1
assay=$2

input_directory="$(pwd)/02-demultiplexed/${assay}"
declare -a SITES=("")
declare -a SAMPLES=({1..5} "WC")

for site in ${SITES[@]}
 do
	mkdir ${input_directory}/${voyage}_${site}
done

for site in ${SITES[@]}
 do
  for sample in ${SAMPLES[@]}
    do
       cp ${input_directory}/${voyage}_${site}_${sample}_${assay}.[12].fq.gz ${input_directory}/${voyage}_${site}
 done
done

mkdir ${input_directory}/Controls
mv ${input_directory}/*EB* ${input_directory}/Controls
mv ${input_directory}/*BC* ${input_directory}/Controls
mv ${input_directory}/NTC* ${input_directory}/Controls
