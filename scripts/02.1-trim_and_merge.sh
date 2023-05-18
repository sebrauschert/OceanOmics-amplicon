#!/bin/bash

# OcOm's 16S: forward_trim = 20, reverse_trim = 22
# OcOm's MiFish: forward_trim = 21, reverse_trim = 27 
# OcOm's COI: forward_trim = 26, reverse_trim = 26 

# Uusage function
usage() {
    echo "Usage: bash scripts/02.1-trim_and_merge.sh -a assay -f forward_trim -r reverse_trim -j max_jobs"
    echo "Options:"
    echo "  -a assay: Specify the assay (16S, MiFish, or COI)"
    echo "  -f forward_trim: Specify the forward trim length (20, 21, or 26)"
    echo "  -r reverse_trim: Specify the reverse trim length (22, 27, or 26)"
    exit 1
}

# Function to trim and concatenate files
process_files() {
    file1=$1
    file2=$2
    suffix=$3
    demux_dir=$4
    forward_trim=$5
    reverse_trim=$6
    prefix=$(basename "$file2" | rev | cut -d '_' -f 2- | rev)

    # Trim reads
    zcat "$file1" | sed "0~2s/^.\{$forward_trim\}//g" > "${demux_dir}/${prefix}_trimmed.${suffix}"
    zcat "$file2" | sed "0~2s/^.\{$reverse_trim\}//g" > "${demux_dir}/${prefix}_trimmed_reverse.${suffix}"

    # Concatenate files
    cat "${demux_dir}/${prefix}_trimmed.${suffix}" "${demux_dir}/${prefix}_trimmed_reverse.${suffix}" | gzip > "${demux_dir}/${prefix}_merged.${suffix}.gz"

    # Remove files that aren't needed any more
    rm "$file1" "$file2" "${demux_dir}/${prefix}_trimmed.${suffix}" "${demux_dir}/${prefix}_trimmed_reverse.${suffix}"
}
export -f process_files

while getopts a:f:r: flag
do
    case "${flag}" in
        a) assay=${OPTARG};;
        f) forward_trim=${OPTARG};;
        r) reverse_trim=${OPTARG};;
        *) usage;;
    esac
done

demux_dir=01-demultiplexed/${assay}
filesFw=(${demux_dir}/*1.fq.gz)
filesRv=(${demux_dir}/*2.fq.gz)

# Loop over read 1 files (there should be two read 1 files for each sample)
for ((i = 0; i < ${#filesFw[@]}; i += 2)); do
    file1=${filesFw[i]}
    file2=${filesFw[i+1]}

    # Trim and concatenate the files
    process_files "$file1" "$file2" "1.fq" "$demux_dir" $forward_trim $reverse_trim &
done
wait

# Loop over read 2 files (there should be two read 2 files for each sample)
for ((i = 0; i < ${#filesRv[@]}; i += 2)); do
    file1=${filesRv[i]}
    file2=${filesRv[i+1]}

    # Trim and concatenate the files
    process_files "$file1" "$file2" "1.fq" "$demux_dir" $forward_trim $reverse_trim &
done
wait

for file in ${demux_dir}/*_merged.*; do 
    mv -- "$file" "${file//_merged/}"; 
done