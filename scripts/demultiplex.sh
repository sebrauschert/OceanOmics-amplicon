#!/usr/bin/env bash

# Demultiplex the reads; As we run different libraries (WGS, metabarcoding and shotgun)
# in one sequencing run, the sampels are demultiplexed based on the library. We need to
# demultiplex the reads by sample before we can continue with DADA2 and blast

# This script used cutadapt v3.2
echo De-multiplexing

input_directory="$(pwd)/00-raw-data"
output_folder="$(pwd)/02-demultiplexed/16S"
read1=$(ls $input_directory/*R1*fastq.gz* | grep -v 'Undetermined*')
read2=$(ls $input_directory/*R2*fastq.gz* | grep -v 'Undetermined*')

# The -g and -G option specify that we are dealing with combinatorial adapters. It further specifies the 5' end, where the adapters are ligated.
# As per cutadapt documentation llumina’s combinatorial dual indexing strategy uses
# a set of indexed adapters on R1 and another one on R2. Unlike unique dual indexes (UDI),
# all combinations of indexes are possible.
# this is another difference: the output will assign the name from the forward and reverse
# reads that were identified with the dual index


echo ${filename}
cutadapt -j 16 \
           -e 0.15 \
           --no-indels \
           -g file:${input_directory}/adapters/RSTM19Feb21_forward_16S.fa  \
           -G file:${input_directory}/adapters/RSTM19Feb21_reverse_16S.fa \
           -o ${output_folder}/{name1}-{name2}.R1.fq.gz \
           -p ${output_folder}/{name1}-{name2}.R2.fq.gz \
           --report=full \
           --minimum-length 1 \
           $read1 $read2


# Need to rename the resulting files to the samples, which I will do as the cutadadapt documentation describes it with mmv
