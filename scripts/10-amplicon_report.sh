#!/bin/bash

voyageID=
assay=
sequencing_run=NA

#..........................................................................................
usage()
{
    printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay; use flag multiple times for multiple assays>\t<string>\n\t\t\t -r <sequencing_run; can be left blank> \n\n";
    exit 1;
}
while getopts v:a:r:w: flag
do
    case "${flag}" in
        v) voyageID=${OPTARG};;
        a) assay+=("$OPTARG");;
        r) sequencing_run=${OPTARG};;
        w) wd=${OPTARG};;
        *) usage;;
    esac
done

if [ "${voyageID}" == ""  ]; then usage; fi
if [ "${assay[1]}" == ""  ]; then usage; fi
if [ "${wd}" == "" ]; then wd=$(pwd); fi

# For the containerised version: if the ANALYSIS path is present,
# change to the ANALYSIS directory
if [ -n "$ANALYSIS" ]
    then cd $ANALYSIS;
fi

# We need to build the assay string in the correct format for the r markdown script (e.g. '16S,MiFish')
assay_rmd_input=
for a in "${assay[@]}"
do
    assay_rmd_input="${assay_rmd_input},${a}"
done
assay_rmd_input="${assay_rmd_input:2}"

# We also need to get the numbers of the samples randomly chosen for creating the quality profile plots
# We can get these numbers from the file names of the plots
random_samples=
sample_plots=($(ls 03-dada2/QC_plots/"${voyageID}"_qualityprofile_Fs_*_"${assay[1]}"_raw.png))
prefix=03-dada2/QC_plots/${voyageID}_qualityprofile_Fs_
suffix=_${assay[1]}_raw.png

for i in "${sample_plots[@]}"
do
    sample="$i"
    sample=${sample#"$prefix"}
    sample=${sample%"$suffix"}
    random_samples="${random_samples},${sample}"
done
random_samples="${random_samples:1}"

Rscript -e "rmarkdown::render('scripts/report/amplicon_report.Rmd',params=list(voyage = '${voyageID}', assays = '${assay_rmd_input}', random_samples = '${random_samples}', sequencing_run = '${sequencing_run}'))"

mv ${wd}/scripts/report/amplicon_report.html ${wd}/06-report