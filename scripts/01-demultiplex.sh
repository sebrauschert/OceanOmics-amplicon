#!/bin/bash
#set -u
#set -o pipefail

# Too avoid too many open files error:
ulimit -S -n 4096

#..........................................................................................

# Take two command line inputs
voyageID=$1
assay=$2

# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate amplicon

# This script used cutadapt v3.2
echo De-multiplexing
input_directory="$(pwd)/00-raw-data"
output_folder="$(pwd)/02-demultiplexed/${assay}"
read1=$(ls $input_directory/*${assay}*R1*fastq.gz* | grep -v 'Undetermined*')
read2=$(ls $input_directory/*${assay}*R2*fastq.gz* | grep -v 'Undetermined*')

# Create the output folder, if it does not already exist
mkdir -p $(pwd)/02-demultiplexed/${assay}

#..........................................................................................
# | The -g and -G option specify that we are dealing with combinatorial adapters.
# | As per cutadapt documentation llumina’s combinatorial dual indexing strategy uses 
# | a set of indexed adapters on R1 and another one on R2. Unlike unique dual indexes (UDI), 
# | all combinations of indexes are possible.
# | this is another difference: the output will assign the name from the forward and reverse
# | reads that were identified with the dual index
#..........................................................................................

echo ${filename}
cutadapt -j 100 \
           -e 0.15 \
           --no-indels \
           -g file:${input_directory}/indices/${voyageID}_${assay}_Fw.fa  \
           -G file:${input_directory}/indices/${voyageID}_${assay}_Rv.fa \
           -o ${output_folder}/{name1}-{name2}.R1.fq.gz \
           -p ${output_folder}/{name1}-{name2}.R2.fq.gz \
           --report=full \
           --minimum-length 1 \
           $read1 $read2

#..........................................................................................
