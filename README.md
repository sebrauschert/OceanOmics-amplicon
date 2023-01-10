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
04-DADA2.R                             - to trim the reads and create an amplicon sequencing variant table; this can be executed in pooled, site specific and fixed error modes
05-run_LULU.sh                         - LULU is a tool to perform post-clustering curation of DNA amplicon data; also creates the final phyloseq object 
06-run_blast.sh                        - to query the NCBI nt and taxa database
07-run_LCA.sh                          - to retrieve the lowest common ancestor for the 16S and MiFish blast results
08-Decontam.R                          - sample decontamination; currently not automated
09-create_phyloseq_object.R            - create the phyloseq R package object for ease of analysis
LCA_filter_nt_only.R                   - filter blast nt results 
Reorganise.sh                          - for site specific analysis: organise samples by site
blast-16S-MiFish.py                    - to query a custom 16S fish database and the MiFish database from here: http://mitofish.aori.u-tokyo.ac.jp/download.html
ecology_plots.R                        - phyloseq based ecology plots for initial alpha and beta diversity 
run_blastnt.sh                         - blast against NCBI nt database

LULU                                   - all functions required for LULU: 01-get_lineage.sh, 02-merge_lineage_with_LCA.R, 03-lulu_create_match_list.sh, 04-LULU.R and 05-create_phyloseq_object.R
dada                                   - three scripts: dada2_pooled.R, dada2_site_error_fixed.R and dada2_site_spec_error.R; called by the 04-DADA2.R code
LCA                                    - LCA scripts and dependencies from eDNAflow pipeline 
```

## Dependencies

### Install `miniconda`

To run this pipeline as smoothly as possible, please install `miniconda` on your system, as per the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) 


### `conda` environments

This repository comes with a `env` folder, which allows to set up three different `conda` environments

- `renv` for a version controlled `R` environment including the `renv` package
- `amplicon` for all utilities required, e.g. `cutadapt`and `seqkit`
- `taxonkit` for taxonomy-related tasks
- `pytaxonkit`dependencies for the python script for the  16S and MiFish blast
-  `blastn` taxonomic annotation of ASVs

To create those environments, first install miniconda, use that to install mamba:

```
conda install -n base mamba
```

Then run the following:

```
mamba env create -f env/renv_environment.yml
mamba env create -f env/amplicon_environment.yml
mamba env create -f env/taxonkit.yml
mamba env create -f env/pytaxonkit.yml
mamba env create -f env/blast-2.12.0.yml
```

There are alternative yml files in the folder for 'general' environments outside of OceanOmics:

```
conda env create -f env/amplicon_environment.general.yml
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

The taxdump-DB is about 500MB in size so make sure you have that space available in your home-directory.

It is possible to extract the data elsewhere and then give taxonkit the path. See the `--data-dir` flag in `computeLCA.py`.

## Databases
Include or script to download? 16S needs Mike Bunce permission.

### NCBI nt
To download the NCBI nt database, we first need to set up the blast conda environment, then download the nt database and the taxdb:

``` shell
conda activate blast-2.12.0
mkdir ncbi-nt
cd ncbi-nt

# nt database
update_blastdb.pl --decompress nt 

# download taxdb in the same folder (ncbi-nt)
wget ftp://ftp.ncbi.nlm.gov/blast/db/taxdb.tar.gz

# Extract
tar xzvf taxdb.tar.gz
```

### 16S
### MiFish

``` shell
mkdir MiFishDB
cd MiFishDB
wget http://mitofish.aori.u-tokyo.ac.jp/files/mitogenomes.zip
unzip mitogenomes.zip
rm *genes.fa
cat *.fa > MiFishDB.fasta
```

## How To

### Bash scripts

#### Set up the analysis environment

Firstly, download this GitHub repository to your local file system, and cd into the repository's folder.

Then to set up the local folder-structure:

```
bash scripts/00-setup.sh myFirstProject
```

This will create a new folder containing subfolders for each step of the pipeline. The folder will be named myFirstProject_Amplicon_YOURUSERNAME in the folder you're currently in.
The folder contains subfolders for all subsequent steps of the pipeline:

```
.
|--- 00-raw-data
|    |--- indices
|--- 01-QC
|--- 02-demultiplexed
|--- 03-dada2
|    |--- QC_plots
|    |--- tmpfiles 
|    |--- errorModel
|--- 04-taxa
|    |--- blast_out
|    |---LCA_out
|--- 05-LULU
|--- 06-report
|--- scripts
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
Rscript scripts/04-DADA2.R --voyage <VOYAGE_ID> \
                           --assay <ASSAY_ID> \
                           --option <pooled, site or fixed>
```

The final result are quality plots before and after read quality trimming, dereplicated reads, merged paired end reads with no chimeras, and .Rdata files for each step in case on step crashes. The end result is an amplicon sequence variant (ASV) table and a fasta of the ASV sequences. 

All results will be in 03-dada2/

#### LULU: post-clustering curation of DNA amplicon data

``` 

```

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
Adam Bennett
