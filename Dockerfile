FROM rocker/tidyverse:4
LABEL maintainer="Seb Rauschert <srauschert@minderoo.org>"
LABEL description="This Dockerfile sets up the Minderoo OceanOmics amplicon pipeline for marine fish eDNA analysis"
MAINTAINER Seb Rauschert <srauschert@minderoo.org>


#=====================================================================================================
# 1. Install essentials, including GNU parallel and mmv, both necessary for the pipeline
# A lot of the below are necessary for the succesful installation of the dada2 package in R
#=====================================================================================================
RUN apt-get update && \
    apt-get install -y build-essential  && \
    apt-get install -y wget && \
    apt-get install -y mmv && \
    apt-get install -y parallel && \
    sed -i 's/# deb-src/deb-src/' /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y liblzma-dev && \
    apt-get -y install libbz2-dev && \
    apt-get -y install libglpk-dev && \
    apt-get -y build-dep libcurl4-gnutls-dev && \
    apt-get -y install libcurl4-gnutls-dev && \
    apt-get -y install datalad && \
    apt-get -y install tree && \
    apt-get -y install pandoc && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get -y autoremove

#=====================================================================================================
# 2. Install conda and set up environments
# Easiest to avoid dependency hell on the system by setting up conda environments.
# With the .yml files this also makes the setup reproducible.
#=====================================================================================================
# from https://fabiorosado.dev/blog/install-conda-in-docker/

# Official minoconda: https://hub.docker.com/r/continuumio/miniconda3/dockerfile
# Install miniconda
ENV CONDA_DIR /opt/conda
ARG CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH
ARG PATH=$CONDA_DIR/bin:$PATH

# Create conda environments
RUN mkdir /opt/conda_env
ADD env /opt/conda_env

RUN conda update -n base -c defaults conda && \
    conda env create -f /opt/conda_env/blast-2.12.0.yml && \
    conda env create -f /opt/conda_env/amplicon_environment.yml && \
    conda env create -f /opt/conda_env/pytaxonkit.yml

# Remove the .yml files
RUN rm -rf /opt/conda_envs


#======================================================================================================
# 3. R RELATED SETUP
# We have a
#======================================================================================================
# install R packages required
RUN mkdir -p /opt/software/setup/R
ADD r_dependencies/install_packages.R /opt/software/setup/R
RUN Rscript /opt/software/setup/R/install_packages.R
RUN rm -rf /opt/software/setup/R/install_packages.R


#======================================================================================================
# 4. Add all pipeline scripts to the container
# With all files added to the container, we do not need to have them on file on the host system.
# Need to copy them from the GitHub repository and make them executable
#======================================================================================================

RUN mkdir /opt/amplicon_pipeline
ADD scripts /opt/amplicon_pipeline

RUN chmod +x /opt/amplicon_pipeline/00-setup.sh && \
    chmod +x /opt/amplicon_pipeline/01-demultiplex.sh && \
    chmod +x /opt/amplicon_pipeline/02-rename_demux.sh && \
    chmod +x /opt/amplicon_pipeline/03-seqkit_stats.sh  && \
    chmod +x /opt/amplicon_pipeline/04.1-docker_DADA2.sh && \
    chmod +x /opt/amplicon_pipeline/05-run_LULU.sh  && \
    chmod +x /opt/amplicon_pipeline/06-run_blast.sh  && \
    chmod +x /opt/amplicon_pipeline/07-run_LCA.sh  && \
    chmod +x /opt/amplicon_pipeline/07.1-LCA_filter_nt_only.R && \
    chmod +x /opt/amplicon_pipeline/08-Decontam.R  && \
    chmod +x /opt/amplicon_pipeline/09-create_phyloseq_object.R  && \
    chmod +x /opt/amplicon_pipeline/10-amplicon_report.sh  && \
    chmod +x /opt/amplicon_pipeline/blast/blast-16S-MiFish.py  && \
    chmod +x /opt/amplicon_pipeline/blast/run_blastnt.sh  && \
    chmod +x /opt/amplicon_pipeline/ecology_plots.R  && \
    chmod +x /opt/amplicon_pipeline/Reorganise.sh  && \
    chmod +x /opt/amplicon_pipeline/LCA/runAssign_collapsedTaxonomy.py && \
    chmod +x /opt/amplicon_pipeline/LCA/working_function.py && \
    chmod +x /opt/amplicon_pipeline/LULU/01-lulu_create_match_list.sh && \
    chmod +x /opt/amplicon_pipeline/LULU/02-LULU.R && \
    chmod +x /opt/amplicon_pipeline/download_nt.sh && \
    chmod +x /opt/amplicon_pipeline/download_mitofish.sh


ENV PATH=/opt/amplicon_pipeline:$PATH
ARG PATH=/opt/amplicon_pipeline:$PATH

# Create a mount directory to write to
RUN mkdir -p /mnt/scratch

#==================================================================================================================
# Making sure that docker starts in the working directory that we mount to the current working directory
# so that we do not need to include an additional path into the scripts. That way, we can add the Dockerfile to the
# pipeline GitHub without modifying the scripts itself, to build the container new every time we update the repo
#
# Therefore, we also set CODE and ANALYSIS environment variables, to make sure the Rscripts are working
#==================================================================================================================

ENV CODE=/opt/amplicon_pipeline/
ENV ANALYSIS=/mnt/scratch/

WORKDIR /mnt/scratch
