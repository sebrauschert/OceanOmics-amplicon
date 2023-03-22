#................................................................
# Installation of required R packages for the amplicon pipeline
#
# Seb Rauschert
# 16/11/2022
#................................................................

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")

BiocManager::install("phyloseq")

BiocManager::install("DECIPHER")

install.packages("phangorn")

BiocManager::install("Biostrings")

install.packages("devtools", repos = 'http://cran.rstudio.com/')

devtools::install_github("tobiasgf/lulu")

install.packages("optparse")

install.packages("tidyverse")

install.packages("png")

install.packages("ggplot2")

install.packages("scales")

install.packages("reshape2")

install.packages("tibble")