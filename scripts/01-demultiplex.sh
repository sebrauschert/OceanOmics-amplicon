#!/bin/bash
#set -u
#set -o pipefail


# USAGE
# bash 01-demultiplex.sh -v <project/voyage ID> \
#                        -a <assay; use teh flag multiple timed for multiple assays, e.g. -a 16S -a MiFish> \
#                        -c <number of cores, default 50>

voyageID=
#assay=
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
#if [ "${assay[@]}" == ""  ]; then usage; fi

# Too avoid too many open files error:
ulimit -S -n 4096


# load amplicon environment
eval "$(conda shell.bash hook)"
conda activate amplicon

set -x
echo 'Writing logs to logs/01-demultiplex.log'
exec 1>logs/01-demultiplex.log 2>&1

for a in "${assay[@]}"
do
    # This script uses cutadapt v4.1
    echo De-multiplexing
    input_directory="$(pwd)/00-raw-data"
    output_folder="$(pwd)/01-demultiplexed/${a}"
    read1=$(ls $input_directory/*${a}*R1*fastq.gz* | grep -v 'Undetermined*')
    read2=$(ls $input_directory/*${a}*R2*fastq.gz* | grep -v 'Undetermined*')

    # Create the output folder, if it does not already exist
    mkdir -p $(pwd)/01-demultiplexed/${a}

    #..........................................................................................
    # | The -g and -G option specify that we are dealing with combinatorial adapters.
    # | As per cutadapt documentation llumina’s combinatorial dual indexing strategy uses
    # | a set of indexed adapters on R1 and another one on R2. Unlike unique dual indexes (UDI),
    # | all combinations of indexes are possible.
    # | this is another difference: the output will assign the name from the forward and reverse
    # | reads that were identified with the dual index
    # |
    # |the '^' in front of file (^file:) means that we anchor the tags to the beginning of the read!
    #..........................................................................................


    cutadapt -j ${cores} \
        -e 0.15 \
        --no-indels \
        -g ^file:${input_directory}/indices/${voyageID}_${a}_Fw.fa  \
        -G ^file:${input_directory}/indices/${voyageID}_${a}_Rv.fa \
        -o ${output_folder}/{name1}-{name2}.R1.fq.gz \
        -p ${output_folder}/{name1}-{name2}.R2.fq.gz \
        --report=full \
        --minimum-length 1 \
        $read1 $read2
done
#..........................................................................................
exit
exitstatus=$?
if [[ $exitstatus != 0 ]]; then
	echo $exitstatus
fi
