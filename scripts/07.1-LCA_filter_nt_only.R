#........................................................................
# ANALYSIS OF AMPLICON DATA: Filtering LCA results for NCBI nt results only
#........................................................................

# Usage: Rscript scripts/LCA_filter.R -v RSV5 -a 16S


suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(RColorBrewer)) 
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))

# this is necessary for the docker version of this script
if(Sys.getenv("ANALYSIS") == ""){

  next

}else{

  setwd(Sys.getenv("ANALYSIS"))

}


# Define options for command line
option_list = list(
  make_option(c("-v", "--voyage"), action="store", default=NA, type='character',
              help="voyage identifier code"),
  make_option(c("-a", "--assay"), action="store", default=NA, type='character',
              help="assay, e.g. '16S' or 'MiFish"))

opt = parse_args(OptionParser(option_list=option_list))

voyage <- opt$voyage
assay  <- opt$assay

# Read in LCA results
lca_results <- read_delim(paste0('05-taxa/LCA_out/',voyage,'_',assay,'_nt_LCA.tsv'))

# View classes for contaminants
print(paste0("Class detected: ", unique(lca_results$class)))

# Check for duplicate ASVs after the LCA
print(paste0("Duplicate ASVs: ", lca_results$OTU[duplicated(lca_results$OTU)]))

# Filter for only chordates (removes all bacterial and euykaryote contaminants, except for mammals)
lca_filtered <- lca_results %>%
  filter(!(domain %in% lca_results$domain[grep(NA, (lca_results$domain))])) %>%
  filter(!(domain %in% lca_results$domain[grep("Bacteria", (lca_results$domain))])) %>%
  filter(!(domain %in% lca_results$domain[grep("Archaea", (lca_results$domain))])) %>%
  filter(!(domain %in% lca_results$domain[grep("dropped", (lca_results$domain))])) %>%
  filter(!(phylum %in% lca_results$phylum[grep(NA, (lca_results$phylum))])) %>%
  filter(!(phylum %in% lca_results$phylum[grep("Chlorophyta", (lca_results$phylum))])) %>%
  filter(!(phylum %in% lca_results$phylum[grep("Bacillariophyta", (lca_results$phylum))])) %>%
  filter(!(phylum %in% lca_results$phylum[grep("dropped", (lca_results$phylum))])) %>%
  filter(!(class %in% lca_results$class[grep(NA, (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Mammalia", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Pelagophyceae", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Cryptophyceae", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Aves", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Caudoviricetes", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Florideophyceae", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Asteroidea", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Ophiuroidea", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("Hexanauplia", (lca_results$class))])) %>%
  filter(!(class %in% lca_results$class[grep("dropped", (lca_results$class))]))

# check classes after filtering - should only be Actinopteri and/or Chondrichthyes
print(paste0("Classes detected after filtering: ", unique(lca_filtered$class)))

# Check for duplicate ASVs (shouldn't be any after decontamination)
print(paste0("Duplicate ASVs after filtering: ", lca_filtered$OTU[duplicated(lca_filtered$OTU)]))

# Write new csv
write.csv(lca_filtered, file = paste0("05-taxa/LCA_out/LCA_filtered_",voyage,"_",assay,"_nt.csv"))
