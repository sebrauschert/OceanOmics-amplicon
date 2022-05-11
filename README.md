# OceanOmics Amplicon Pipeline


<p align="center">
  <img width="330" height="300" src="img/OceanOmics.png">
</p>


## Overview
This repository contains the amplicon sequencing pipeline for taxonomic annotation of the samples collected in the OceanOmics project. It contains the `dada2` pipeline and 
a `blastn` query for taxonomic identification.

The scripts contained herein are:

```
01. 00-setup.sh                            - to set up the analysis folder structure
02. 01-demultiplex.sh                      - to demultiplex the amplicon data added to the 00-raw-data folder
03. 02-rename_demux.sh                     - to rename the demultiplexed reads to their respective sample names; a index to sample name mapping file needs to be created
04. 03-seqkit_stats.sh                     - to create read statistics for QC checks of the demultiplexed reads
05. 04-DADA2.R                             - to trim the reads and create an amplicon sequencing variant table
06. 05-blastn.sh                           - to query the NCBI nt and taxa database
07. 06-LCA/runAssign_collapsedTaxonomy.py  - a custom script from the eDNAflow pipeline for lowest common ancestor analysis of the taxonomically annotated ASVs
08. 07-ecology_plots.R                     - phyloseq based ecology plots for initial alpha and beta diversity 
```

## Dependencies

### Install `miniconda`

To run this pipeline as smoothly as possible, please install `miniconda` on your system, as per the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) 


### `conda` environments
This repository comes with a `env` folder, which allows to set up three different `conda` environments

- `datalad` environment, to track data and analysis
- `renv` for a version controlled `R` environment including the `renv` package
- `ampplicon` for all utilities required, e.g. `cutadapt`and `seqkit`

To create those environments, first install miniconda end then run the following:

```
conda env create -f env/datalad_environment.yml
conda env create -f env/renv_environment.yml
conda env create -f env/amplicon_environment.yml
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

## How To

### Bash scripts

#### Set up the analysis environment

Firstly, download this GitHub repository to your local file system, and cd into the repository's folder.

Then to set up the local datalad project:

```
bash scripts/00-setup.sh myFirstProject
```

This will create a new folder containing a datalad project. It will be named myFirstProject_Amplicon_YOURUSERNAME in the folder you're currently in.
The folder contains subfolders for all subsequent steps of the pipeline.

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

### `Nextflow`

## Authors and contributors
Jessica Pearce  
Sebastian Rauschert  
Priscila Goncalves  
