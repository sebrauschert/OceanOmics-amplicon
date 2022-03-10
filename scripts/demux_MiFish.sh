#!/bin/bash
mkdir -p 02-demultiplexed

# This script used cutadapt v3.2
echo De-multiplexing
input_directory="$(pwd)/00-raw-data"
output_folder="$(pwd)/01-demultiplexed/MiFish"
read1=$(ls $input_directory/*R1*fastq.gz* | grep -v 'Undetermined*')
read2=$(ls $input_directory/*R2*fastq.gz* | grep -v 'Undetermined*')

# The -g and -G option specify that we are dealing with combinatorial adapters.
# As per cutadapt documentation llumina’s combinatorial dual indexing strategy uses 
# a set of indexed adapters on R1 and another one on R2. Unlike unique dual indexes (UDI), 
# all combinations of indexes are possible.
# this is another difference: the output will assign the name from the forward and reverse
# reads that were identified with the dual index
#
# Change log: I changed the -a and -A to -g and -G; This seems to pick up signals and I will try it afterwards with the -a -A option again
# It seems, however, that the -g and -G is correct. Need to clarify this
echo ${filename}
cutadapt -j 16 \
           -e 0.15 \
           --no-indels \
           -g file:${input_directory}/adapters/RSTM19Feb21_forward_MiFish.fa  \
           -G file:${input_directory}/adapters/RSTM19Feb21_reverse_MiFish.fa \
           -o ${output_folder}/{name1}-{name2}.R1.fq.gz \
           -p ${output_folder}/{name1}-{name2}.R2.fq.gz \
           --report=full \
           --minimum-length 1 \
           $read1 $read2

# Create a stats table of reads and read length per sample
seqkit stats -T *fq.gz > demux_seqkit_stats_MiFish.txt


# Need to rename the resulting files to the samples, which I will do as the cutadadapt documentation describes it with mmv
