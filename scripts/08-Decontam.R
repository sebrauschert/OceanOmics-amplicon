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
              help="assay, e.g. '16S' or 'MiFish"),
  make_option(c("-o", "--option"), action="store", default=NA, type='character',
              help="nt or custom blast database"))  

opt = parse_args(OptionParser(option_list=option_list))

voyage <- opt$voyage
assay  <- opt$assay
option <- opt$option

# get vector of control files that end in .1.fq.gz
suffix <- paste0("_", assay, ".1.fq.gz")
controls <- list.files(paste0("./01-demultiplexed/", assay, "/Controls/"), pattern = paste0("*", suffix))

# the water controls might be in the site folders
water_suffix <- paste0("WC", suffix)
water_controls <- list.files(paste0("./01-demultiplexed/", assay, "/"), pattern = paste0("*", water_suffix), recursive = TRUE)                       
water_controls <- basename(water_controls)

# concatenate the two vectors and remove the suffix
controls <- c(controls, water_controls)
controls <- sub(suffix, "", controls)

#......................................................................................
# WE CALL THE SCRIPT WE NEED BASED ON THE OPTIONS INPUT
#......................................................................................

# Run the analysis by executing the function above

 if(option == "nt"){
   
   # Read in filtered LCA results
   lca_tab <- read_csv(paste0("05-taxa/LCA_out/LCA_filtered_", voyage, "_", assay, "_nt.csv"))
   
   # Mark all potential contaminant ASV sequences in new column
   lca_tab$Contam <- "False"
   
   # All ASV sequences identified in any control samples
   for (i in controls){
     lca_tab$Contam[lca_tab[i] >0] <- "True"
   }
   
   lca_tab <- lca_tab %>%
     relocate(OTU, .before = domain) %>%
     rename(ASV = OTU)
   
   # Write final output with contam labels
   write_csv(lca_tab, file = paste0("05-taxa/", voyage, "_", assay, "_contam_table_nt.csv"))
   
   # Create final file with no contaminants
   nocontam <- read_csv(file = paste0("05-taxa/", voyage, "_", assay, "_contam_table_nt.csv"))
   nocontam$Contam
   nocontam <- subset(nocontam, Contam=="FALSE")
   
   nocontam <- nocontam %>% 
     select(where(~ any(. != 0)))
   
   write_csv(nocontam, file = paste0("05-taxa/", voyage, "_", assay, "_nocontam_nt.csv"))
  
   
   
 }

 if(option == "custom"){
   
   # Read in filtered LCA results
   lca_tab <- read_delim(paste0("05-taxa/LCA_out/", voyage, "_", assay, "_LCA.tsv"))
   
   # Mark all potential contaminant ASV sequences in new column
   lca_tab$Contam <- "False"
   
   # All ASV sequences identified in any control samples
   for (i in controls){
     lca_tab$Contam[lca_tab[i] >0] <- "True"
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
