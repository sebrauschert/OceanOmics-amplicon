# OceanOmics Amplicon Pipeline


<p align="center">
  <img width="330" height="300" src="img/OceanOmics.png">
</p>


## Overview
This repository contains the amplicon sequencing pipeline for taxonomic annotation of the samples collected in the OceanOmics project. It contains the `dada2` pipeline and 
a `blastn` query for taxonomic identification.

The scripts contained herein are:

```
00-setup.sh                            - to set up the analysis folder structure
01-demultiplex.sh                      - to demultiplex the amplicon data added to the 00-raw-data folder
02-rename_demux.sh                     - to rename the demultiplexed reads to their respective sample names; a index to sample name mapping file needs to be created
03-seqkit_stats.sh                     - to create read statistics for QC checks of the demultiplexed reads
04-DADA2.R                             - to trim the reads and create an amplicon sequencing variant table
05-blastn.sh                           - to query the NCBI nt and taxa database
06-blast-16S-MiFish.py                 - to query a custom 16S fish database and the MiFish database from here: http://mitofish.aori.u-tokyo.ac.jp/download.html
07-custom-lca.py                       - to retrieve the lowest common ancestor for the 16S and MiFish blast results
08-LCA/runAssign_collapsedTaxonomy.py  - a custom script from the eDNAflow pipeline for lowest common ancestor analysis of the taxonomically annotated ASVs
09-ecology_plots.R                     - phyloseq based ecology plots for initial alpha and beta diversity 
```

## Dependencies

### Install `miniconda`

To run this pipeline as smoothly as possible, please install `miniconda` on your system, as per the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) 


### `conda` environments

This repository comes with a `env` folder, which allows to set up three different `conda` environments

- `datalad` environment, to track data and analysis
- `renv` for a version controlled `R` environment including the `renv` package
- `amplicon` for all utilities required, e.g. `cutadapt`and `seqkit`
- `taxonkit` for taxonomy-related tasks
- `pytaxonkit`dependencies for the python script for the  16S and MiFish blast

To create those environments, first install miniconda end then run the following:

```
conda env create -f env/datalad_environment.yml
conda env create -f env/renv_environment.yml
conda env create -f env/amplicon_environment.yml
conda env create -f env/taxonkit.yml
conda env create -f env/pytaxonkit.yml
```

There are alternative yml files in the folder for 'general' environments outside of OceanOmics:

```
conda env create -f env/datalad_environment.general.yml
```

### Install `mmv`

This pipeline requires the linux utility `mmv` installed on your system.

Please install it with:

```
sudo apt-get update
sudo apt-get install mmv
```

### Install `taxonkit`

Following the above conda instructions, we now have a conda environment called `taxonkit`. taxonkit expects the NCBI taxdump in ~/.taxonkit:
TODO: potentially removing this dependency and use pytaxonkit exclusively?
```
mkdir ~/.taxonkit
cd ~/.taxonkit
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzvf taxdump.tar.gz
```
It's about 500MB in size.

It is possible to extract the data elsewhere and then give taxonkit the path. See the `--data-dir` flag in `computeLCA.py`.

## Databases
Include or script to download? 16S needs Mike Bunce permission.
### 16S
### MiFish

## How To

### Bash scripts

#### Set up the analysis environment

Firstly, download this GitHub repository to your local file system, and cd into the repository's folder.

Then to set up the local datalad project:

```
conda activate datalad
bash scripts/00-setup.sh myFirstProject
```

This will create a new folder containing a datalad project. It will be named myFirstProject_Amplicon_YOURUSERNAME in the folder you're currently in.
The folder contains subfolders for all subsequent steps of the pipeline:

```

.
├── 00-raw-data
│   └── indices
├── 01-QC
├── 02-demultiplexed
│   └── sample_names
├── 03-dada2
│   └── QC_plots
├── 04-taxa
├── 05-report
└── scripts
    └── LCA
```

The pipeline expects fastq of your amplicon sequencing data in the 00-raw-data/ folder. The files should be named `*R1*fastq.gz` and `*R2*.fastq.gz`.

#### (optional) Demultiplexing

Some amplicon data needs to be demultiplexed. We store indices for demultiplexing in a voyageID/assayID structure. For each voyage and each assay we need two files in `00-raw-data/indices`: `voyageID_assayID_Fw.fa` and `voyageID_assayID_Rv.fa`. Please add these manually.

Then, to run the demultiplexing with those indices:
```
conda activate amplicon
bash scripts/01-demultiplex.sh Voyage1 assay1
```

This will run cutadapt with all fastq files in 00-raw-data and the indices for this particular voyage/assay combination. The demultiplexed reads will be in 02-demultiplexed/voyageID/assayID


Now we need to rename the output files to correspond to their sample IDs.

```
conda activate amplicon
bash scripts/02-rename_demux.sh
```

For each voyage's assay this script will rename the demultiplexed fastq files according to their sample ID.

You can calculate statistics for each assay using the seqkit_stats script:

```
conda activate amplicon
bash scripts/03-seqkit_stats.sh
```

This will generate a txt file of QC statistics for each assay and voyage in the folder 01-QC/


#### Amplicon sequence variants via DADA2 


We run the R-script DADA2 to calculate amplicon sequence variants (ASVs).

```
conda activate renv
Rscript scripts/04-DADA2.R
```
This scripts has some hardcoded variable names that need to be changed - notably `voyages` and `assays` at the bottom of the script.

This will run DADA2 for every separately for each voyage and each assay. The final result are quality plots before and after read quality trimming, dereplicated reads, merged paired end reads with no chimeras, and .Rdata files for each step in case on step crashes. The end result is an amplicon sequence variant (ASV) table and a fasta of the ASV sequences. 

All results will be in 03-dada2/

#### Taxonomic assignment via blastn

This step uses blastn to find taxonomic hits for the ASV fasta sequences.

This needs to be run for each voyage and each assay, here an example requesting 12 CPUs and writing to the output file 04-taxa/tabular_output_voyage1assay1.tsv

```
conda activate amplicon
bash scripts/05-blastn.sh 03-dada2/voyage1_assay1.fa 12 voyage1assay1
```

#### Taxonomic assignment via blastn: 16S and MiFish database
For this we need to activate the `pytaxonkit` environment and execute the script `06-blast-16S-MiFish.py`.

``` 
conda activate pytaxonkit
bash scripts/06-blast-16S-MiFish.py \
          --dada2_file [Path to the dada2 fasta ASV sequence file] \
          --out_path [path to the folder where the output shall be saved; default current working directory] \
          --database [either 16S or MiFish]
```

#### Finding the last common ancestor (LCA) for each query ASV

Now we try to find the 'best' hit for each query by merging all hits into their LCA species.

First, we filter the blast tabular output:

```
awk '{if ($7 > 90) print}' all_results.tsv > all_results.90perc.tsv
grep -v -e uncultured -e Uncultured -e chloroplast -e Unidentified -e unidentified all_results.90perc.tsv > all_results.90perc.noUnculturedUnidentifiedChloroplast.tsv
```

Then, to make a table of the lineage and hit for each query's LCA:

```
conda activate taxonkit
python computeLCA.py all_results.90perc.noUnculturedUnidentifiedChloroplast.tsv > all_results.90perc.noUnculturedUnidentifiedChloroplast.LCAs.tsv
```

Then we can search for fish in these LCAs.

```
python findFish.py all_results.90perc.noUnculturedUnidentifiedChloroplast.LCAs.tsv > all_results.90perc.noUnculturedUnidentifiedChloroplast.LCAs.FishOnly.tsv
# give me the number of hits
wc -l all_results.90perc.noUnculturedUnidentifiedChloroplast.LCAs.FishOnly.tsv
```

### `Nextflow`

## Authors and contributors
Jessica Pearce  
Sebastian Rauschert  
Priscila Goncalves  
Philipp Bayer
