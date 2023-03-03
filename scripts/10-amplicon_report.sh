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
while getopts v:a:r: flag
do
    case "${flag}" in
        v) voyageID=${OPTARG};;
        a) assay+=("$OPTARG");;
        r) sequencing_run=${OPTARG};;
        *) usage;;
    esac
done

if [ "${voyageID}" == ""  ]; then usage; fi
if [ "${assay[1]}" == ""  ]; then usage; fi

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

# For the containerised version: if the CODE path is present,
# change to the CODE directory
if [ -n "$CODE" ]
    then cd $CODE;
fi

Rscript -e "rmarkdown::render('report/amplicon_report.Rmd',params=list(voyage = '${voyageID}', assays = '${assay_rmd_input}', random_samples = '${random_samples}', sequencing_run = '${sequencing_run}'))"

mv report/amplicon_report.html ../06-report
