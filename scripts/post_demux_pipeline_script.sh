usage()
{
    printf "Usage: $0 -v <voyageID>\t<string>\n\t\t\t -a <assay>\t<string>\n\t\t\t -d <demultiplex pipeline output directory>\t<string>\n\n";
    exit 1;
}
while getopts v:a:d: flag
do
    case "${flag}" in
        v) voyageID=${OPTARG};;
        a) assay=${OPTARG};;
        d) demux_out_dir=${OPTARG};;
        *) usage;;
    esac
done
if [ "${voyageID}" == ""  ]; then usage; fi

mkdir -p 00-raw-data/indices
mkdir -p 01-demultiplexed/${assay}
mkdir 02-QC
mkdir 06-report

cp ${demux_out_dir}/renamed_fqs/missing_samples.txt 02-QC/missing_samples_${assay}.txt
cp ${demux_out_dir}/seqkit_stats/assigned_seqkit_stats.txt 02-QC/Sample_statistics_${voyageID}_${assay}_sample_level.txt
cp ${demux_out_dir}/seqkit_stats/unknown_seqkit_stats.txt 02-QC/Sample_statistics_${voyageID}_${assay}_sample_level_unknown.txt
cp ${demux_out_dir}/seqkit_stats/unnamed_seqkit_stats.txt 02-QC/Sample_statistics_${voyageID}_${assay}_sample_level_unnamed.txt
cp ${demux_out_dir}/concat_fqs/* 01-demultiplexed/${assay}
cp ${demux_out_dir}/valid_input/valid_metadata.csv 06-report/${voyageID}_metadata.csv
cp ${demux_out_dir}/index_file/indices.csv 00-raw-data/indices/${voyageID}_${assay}_indices.csv
cp ${demux_out_dir}/demux_dependencies/${assay}_fw.fa 00-raw-data/indices/${voyageID}_${assay}_fw.fa
cp ${demux_out_dir}/demux_dependencies/${assay}_rv.fa 00-raw-data/indices/${voyageID}_${assay}_rv.fa
cp ${demux_out_dir}/demux_dependencies/${assay}_sample_name_rename_pattern.txt 00-raw-data/indices/${voyageID}_${assay}_sample_name_rename_pattern.txt

echo post demux pipeline script finished