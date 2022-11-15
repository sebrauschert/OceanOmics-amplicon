#===========================================================================
# Create phyloseq object (with phyloseq)
# based on here: https://joey711.github.io/phyloseq/import-data.html
# Seb Rauschert
# 22/07/2022
# Modified: Jessica Pearce
# 04/11/2022

# Usage: Rscript scripts/08-create_phyloseq_object.R -v RSV5 -a 16S -o custom
#===========================================================================

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(RColorBrewer)) 
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

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


# Run the analysis by executing the function above
if(option == "nt"){
  
  contam_tab <- read_csv(paste0("05-taxa/",voyage,"_",assay,"_contam_table_nt.csv"))
  dada_tab <- read_tsv(paste0("03-dada2/",voyage,"_final_table_",assay,".tsv"))
  
  asv_seq <- dada_tab %>%
    select(ASV, ASV_sequence)
  
  final_tab <- merge(contam_tab, asv_seq, by = "ASV", all.x = TRUE)
  write_csv(final_tab, paste0("05-taxa/",voyage,"_",assay,"_final_tab_nt.csv"))
  
# Specify paths to data 
otu_mat    <- paste0('05-taxa/',voyage,'_',assay,'_final_tab_nt.csv')
tax_mat    <- paste0('05-taxa/',voyage,'_',assay,'_final_tab_nt.csv')
samples_df <- paste0("06-report/",voyage,"_metadata.csv")

#===========================================================================
# CUSTOM WRANGLING

meta <- read_csv(samples_df)
taxa <- read_csv(tax_mat)
otu  <- read_csv(otu_mat)

# Prepare the metadata
meta             <- as.data.frame(meta)
rownames(meta)   <- meta$`Sample ID`
meta$`Sample ID` <- NULL

# Prepare the taxa data
taxa %>%
  select(ASV, domain, phylum, class, order, family, genus, species, Contam, ASV_sequence) -> taxa
taxa           <- as.data.frame(taxa)
rownames(taxa) <- taxa$ASV
taxa$ASV       <- NULL
taxa           <- as.matrix(taxa)

# Prepare otu data
otu           <- as.data.frame(otu)
rownames(otu) <- otu$ASV
otu[,c('ASV', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'Contam', 'ASV_sequence')] <- list(NULL)

# Create the phyloseq object
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
META = sample_data(meta)

physeq = phyloseq(OTU, TAX, META)
saveRDS(physeq, file = paste0('06-report/',voyage,'_',assay,'_phyloseq_nt.rds'))

#===========================================================================

}

if(option == "custom"){
  
  contam_tab <- read_csv(paste0("05-taxa/",voyage,"_",assay,"_contam_table.csv"))
  dada_tab <- read_tsv(paste0("03-dada2/",voyage,"_final_table_",assay,".tsv"))
  
  asv_seq <- dada_tab %>%
    select(ASV, ASV_sequence)
  
  final_tab <- merge(contam_tab, asv_seq, by = "ASV", all.x = TRUE)
  write_csv(final_tab, paste0("05-taxa/",voyage,"_",assay,"_final_tab.csv"))
  
  
  # Specify paths to data 
  otu_mat    <- paste0('05-taxa/',voyage,'_',assay,'_final_tab.csv')
  tax_mat    <- paste0('05-taxa/',voyage,'_',assay,'_final_tab.csv')
  samples_df <- paste0("06-report/",voyage,"_metadata.csv")
  
  #===========================================================================
  # CUSTOM WRANGLING
  
  meta <- read_csv(samples_df)
  taxa <- read_csv(tax_mat)
  otu  <- read_csv(otu_mat)
  
  # Prepare the metadata
  meta             <- as.data.frame(meta)
  rownames(meta)   <- meta$`Sample ID`
  meta$`Sample ID` <- NULL
  
  # Prepare the taxa data
  taxa %>%
    select(ASV, domain, phylum, class, order, family, genus, species, Contam, ASV_sequence) -> taxa
  taxa           <- as.data.frame(taxa)
  rownames(taxa) <- taxa$ASV
  taxa$ASV       <- NULL
  taxa           <- as.matrix(taxa)
  
  # Prepare otu data
  otu           <- as.data.frame(otu)
  rownames(otu) <- otu$ASV
  otu[,c('ASV', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'Contam', 'ASV_sequence')] <- list(NULL)
  
  # Create the phyloseq object
  OTU = otu_table(otu, taxa_are_rows = TRUE)
  TAX = tax_table(taxa)
  META = sample_data(meta)
  
  physeq = phyloseq(OTU, TAX, META)
  saveRDS(physeq, file = paste0('06-report/',voyage,'_',assay,'_phyloseq.rds'))
  
  #===========================================================================
  
}

### Done!
print(paste0(voyage, " ", assay, " ", option, " phyloseq object successfully created!"))



