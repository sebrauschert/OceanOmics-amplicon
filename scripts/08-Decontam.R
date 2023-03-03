#........................................................................
# ANALYSIS OF AMPLICON DATA: Decontamination of LCA results
#........................................................................

# Currently this script cannot be automated and is voyage specific
# Usage: Rscript scripts/Decontam.R -v RSV5 -a 16S -o custom

suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(RColorBrewer)) 
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

# this is necessary for the docker version of this script
if(Sys.getenv("ANALYSIS") != ""){

  setwd(Sys.getenv("ANALYSIS"))

}

# Define options for command line
option_list = list(
  make_option(c("-v", "--voyage"), action="store", default=NA, type='character',
              help="voyage identifier code"),
  make_option(c("-a", "--assay"), action="store", default=NA, type='character',
              help="assay, e.g. '16S' or 'MiFish"),
  make_option(c("-o", "--option"), action="store", default=NA, type='character',
              help="ocom, nt, or custom blast database"),
  make_option(c("-c", "--control_grep_patterns"), action="store", default=NA, type='character',
              help="Comma separated list of control patterns; e.g., 'WC,EB'"))  

opt = parse_args(OptionParser(option_list=option_list))

voyage                <- opt$voyage
assay                 <- opt$assay
option                <- opt$option
control_grep_patterns <- opt$control_grep_patterns

control_grep_patterns <- str_split_1(control_grep_patterns, ",")

# get vector of demultiplexed files that end in .1.fq.gz then remove the suffux
suffix <- paste0("_", assay, ".1.fq.gz")
demux_files <- list.files(paste0("./01-demultiplexed/", assay, "/"), pattern = paste0("*", suffix))
demux_files <- sub(suffix, "", controls)

#......................................................................................
# WE CALL THE SCRIPT WE NEED BASED ON THE OPTIONS INPUT
#......................................................................................

# Run the analysis by executing the function above


if(option %in% c("nt", "ocom")){
  # we could have nt output, or OcOm database output
  suffix <- option
   
  if (suffix == "nt") {
    # Read in filtered LCA results
    lca_tab <- read_csv(paste0("05-taxa/LCA_out/LCA_filtered_", voyage, "_", assay, "_", suffix, ".csv"))
  } else {
    # Read in filtered LCA results
    lca_tab <- read_tsv(paste0("05-taxa/LCA_out/", voyage, "_", assay, "_", suffix, "_LCA.tsv"))
  }


  # Mark all potential contaminant ASV sequences in new column
  lca_tab$Contam <- "False"
   
  # Flag all ASV sequences identified in control samples
  for (i in demux_files){
    for (pattern in control_grep_patterns) {
      if (grepl(pattern, i, fixed = TRUE)) {
        lca_tab$Contam[lca_tab[i] >0] <- "True"
      }
    }
  }
   
  lca_tab <- lca_tab %>%
    relocate(OTU, .before = domain) %>%
    rename(ASV = OTU)
   
  # Write final output with contam labels
  write_csv(lca_tab, file = paste0("05-taxa/", voyage, "_", assay, "_contam_table_", suffix, ".csv"))
   
  # Create final file with no contaminants
  nocontam <- read_csv(file = paste0("05-taxa/", voyage, "_", assay, "_contam_table_", suffix, ".csv"))
  nocontam$Contam
  nocontam <- subset(nocontam, Contam=="FALSE")

  nocontam <- nocontam %>% 
    select(where(~ any(. != 0)))
   
  write_csv(nocontam, file = paste0("05-taxa/", voyage, "_", assay, "_nocontam_", suffix, ".csv"))   
}



if(option == "custom"){
   
  # Read in filtered LCA results
  lca_tab <- read_delim(paste0("05-taxa/LCA_out/", voyage, "_", assay, "_LCA.tsv"))
   
  # Mark all potential contaminant ASV sequences in new column
  lca_tab$Contam <- "False"
   
  # Flag all ASV sequences identified in control samples
  for (i in demux_files){
    for (pattern in control_grep_patterns) {
      if (grepl(pattern, i, fixed = TRUE)) {
        lca_tab$Contam[lca_tab[i] >0] <- "True"
      }
    }
  }
   
  lca_tab <- lca_tab %>%
    relocate(OTU, .before = domain) %>%
    rename(ASV = OTU)
   
  # Write final output with contam labels
  write_csv(lca_tab, file = paste0("05-taxa/", voyage, "_", assay, "_contam_table.csv"))
   
  # Create final file with no contaminants
  nocontam <- read_csv(file = paste0("05-taxa/", voyage, "_", assay, "_contam_table.csv"))
  nocontam$Contam
  nocontam <- subset(nocontam, Contam=="FALSE")
   
  nocontam <- nocontam %>% 
    select(where(~ any(. != 0)))
   
  write_csv(nocontam, file = paste0("05-taxa/", voyage, "_", assay, "_nocontam.csv"))
}
 
### Done!
print(paste0(voyage, " ", assay, " decontamination done!"))