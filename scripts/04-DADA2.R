#........................................................................
# ANALYSIS OF AMPLICON DATA: DADA2
#........................................................................

suppressPackageStartupMessages(library(dada2)) 
suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(RColorBrewer)) 
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))

# this is necessary for the docker version of this script
if(Sys.getenv("CODE") == ""){

  next

}else{

  setwd(Sys.getenv("CODE"))

}


# All functions for dada2 are in a seperate script now. 
# Seperate functions for the pooled, the site specific and 
# the fixed error rate dada2 analysis

#source("scripts/functions/dada2_functions.R")

# Define options for command line
option_list = list(
  make_option(c("-v", "--voyage"), action="store", default=NA, type='character',
              help="voyage identifier code"),
  make_option(c("-a", "--assay"), action="store", default=NA, type='character',
              help="assay, e.g. '16S' or 'MiFish"),
  make_option(c("-o", "--option"), action="store", default=NA, type='character',
              help="pooled, site or fixed error"),
  make_option(c("-c", "--cores"), action="store", default=20, type='integer',
              help="number of cores, default 20"))

opt = parse_args(OptionParser(option_list=option_list))

voyage <- opt$voyage
assay  <- opt$assay
option <- opt$option
cores  <- opt$cores


#......................................................................................
# WE CALL THE SCRIPT WE NEED BASED ON THE OPTIONS INPUT
#......................................................................................

# Run the analysis by executing the function above
if(option == "pooled"){
  
  source("dada/dada2_pooled.R")
  
}

if(option == "site"){
  
  source("./scripts/dada/dada2_site_spec_error.R")
  
}

if(option == "fixed"){
  
  source("./scripts/dada/dada2_site_error_fixed.R")
  
  
}
    

### Done!
print(paste0(voyage, " ", assay, " DADA2 run finished!"))
