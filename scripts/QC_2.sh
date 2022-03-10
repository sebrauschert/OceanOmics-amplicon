#!/bin/bash


source /data/software/miniconda3/etc/profile.d/conda.sh

#"$(conda shell.bash hook)"
conda activate amplicon

input_directory=$1
assay1=$2
assay2=$3
output_directory=04-qc_2/test/${assay1}
reads=$(ls $input_directory/${assay1}/*f*)

# FastQC to check the quality of the demultiplexed samples;
fastqc $reads \
       -o $output_directory #\
#       -t 16

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

wait # fastqc results need to be available before continuing

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

input_directory=$1
assay1=$2
assay2=$3
output_directory=04-qc_2/test/${assay2}
reads=$(ls $input_directory/${assay2}/*f*)

# FastQC to check the quality of the demultiplexed samples;
fastqc $reads \
       -o $output_directory #\
#       -t 16

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

wait # fastqc results need to be available before continuing

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# MultiQC to create nice output of first assay data

source /data/software/miniconda3/etc/profile.d/conda.sh

#"$(conda shell.bash hook)"
conda activate multiqc

input_directory=04-qc_2/test/${assay1}
output_directory=04-qc_2/test/${assay1}/MultiQC

# MuliQC to check the quality overall samples;
multiqc $input_directory \
       -o $output_directory
       
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

wait # fastqc results need to be available before continuing

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# MultiQC on second assay data

input_directory=04-qc_2/test/${assay2}
output_directory=04-qc_2/test/${assay2}/MultiQC

# MuliQC to check the quality overall samples;
multiqc $input_directory \
       -o $output_directory
