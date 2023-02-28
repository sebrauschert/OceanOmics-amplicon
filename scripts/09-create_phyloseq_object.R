#===========================================================================
# Create phyloseq object (with phyloseq)
# based on here: https://joey711.github.io/phyloseq/import-data.html
# Seb Rauschert
# 22/07/2022
# Modified: Jessica Pearce, Adam Bennett
#
# Usage: Rscript scripts/08-create_phyloseq_object.R -v RSV5 -a 16S -o custom
#===========================================================================

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(Biostrings))

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
  make_option(c("-c", "--cores"), action="store", default=NA, type='integer',
              help="number of cores"),
  make_option(c("-t", "--optimise_tree"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE; setting this to TRUE will optimise the phylogenetic tree and cause the script take longer to run"))


opt = parse_args(OptionParser(option_list=option_list))

voyage   <- opt$voyage
assay    <- opt$assay
option   <- opt$option
cores    <- opt$cores
optimise <- opt$optimise


# Run the analysis by executing the function above
if(option %in% c("ocom", "nt")){
  suffix <- option
  contam_tab <- read_csv(paste0("05-taxa/",voyage,"_",assay,"_contam_table_", suffix, ".csv"))
  dada_tab <- read_tsv(paste0("03-dada2/",voyage,"_final_table_",assay,".tsv"))

  asv_seq <- dada_tab %>%
    select(ASV, ASV_sequence)

  final_tab <- merge(contam_tab, asv_seq, by = "ASV", all.x = TRUE)
  write_csv(final_tab, paste0("05-taxa/",voyage,"_",assay,"_final_tab_", suffix, ".csv"))

# Specify paths to data
otu_mat    <- paste0('05-taxa/',voyage,'_',assay,'_final_tab_', suffix, '.csv')
tax_mat    <- paste0('05-taxa/',voyage,'_',assay,'_final_tab_', suffix, '.csv')
samples_df <- paste0("06-report/",voyage,"_metadata.csv")
tree_mat   <- paste0('05-taxa/',voyage,'_',assay,'_final_tab_', suffix, '.csv')


  #===========================================================================
  # CUSTOM WRANGLING

  meta <- read_csv(samples_df)
  taxa <- read_csv(tax_mat)
  otu  <- read_csv(otu_mat)
  seq_tab  <- read_csv(tree_mat)

  # Prepare the metadata
  meta             <- as.data.frame(meta)

  if (! 'Sample ID' %in% colnames(meta)) {
    stop(paste0("Please make sure ", samples_df, " contains a 'Sample ID' column"));
  }

  rownames(meta)   <- meta$`Sample ID`
  meta$`Sample ID` <- NULL

  # Prepare the taxa data
  taxa["LCA"] <- ""
  taxa %>%
    select(ASV, domain, phylum, class, order, family, genus, species, LCA, Contam, ASV_sequence) -> taxa
  taxa           <- as.data.frame(taxa)

  # We need to fill the LCA column starting with 'species', then 'genus', etc
  # until we reach a taxonomy level that isn't 'na' or 'dropped'
  levels <- c("species", "genus", "family", "order", "class", "phylum", "domain")
  for (row in 1:nrow(taxa)) {
    LCA <- NA
    level_ID <- 1
    while(LCA == "dropped" | is.na(LCA)) {
      LCA <- taxa[row, levels[level_ID]]
      level_ID <- level_ID + 1
    }
    taxa[row, "LCA"] <- LCA
  }

  rownames(taxa) <- taxa$ASV
  taxa$ASV       <- NULL
  taxa           <- as.matrix(taxa)

  # Prepare otu data
  otu           <- as.data.frame(otu)
  rownames(otu) <- otu$ASV
  otu[,c('ASV', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'Contam', 'ASV_sequence')] <- list(NULL)

  # Prepare the phylo tree data
  DNA_set = DNAStringSet(seq_tab$ASV_sequence)
  names(DNA_set) = paste0(seq_tab$ASV)

  alignment = AlignSeqs(DNA_set, anchor=NA, processors=cores)

  phang_align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang_align)
  treeNJ <- NJ(dm)

  fit = pml(treeNJ, data=phang_align)
  fitGTR <- update(fit, k=4, inv=0.2)

  if (optimise) {
    # Use 'optim.pml()' to optimise the model parameters
    fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
  }

  # Create the phyloseq object
  OTU = otu_table(otu, taxa_are_rows = TRUE)
  TAX = tax_table(taxa)
  META = sample_data(meta)
  TREE = phy_tree(fitGTR$tree)

  physeq = phyloseq(OTU, TAX, META, TREE)
saveRDS(physeq, file = paste0('06-report/',voyage,'_',assay,'_phyloseq_', suffix, '.rds'))

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
  tree_mat   <- paste0('05-taxa/',voyage,'_',assay,'_final_tab_', suffix, '.csv')

  #===========================================================================
  # CUSTOM WRANGLING

  meta <- read_csv(samples_df)
  taxa <- read_csv(tax_mat)
  otu  <- read_csv(otu_mat)

  # Prepare the metadata
  meta             <- as.data.frame(meta)

  if (! 'Sample ID' %in% colnames(meta))
  {
    stop(paste0("Please make sure ", samples_df, " contains a 'Sample ID' column"));
  }

  rownames(meta)   <- meta$`Sample ID`
  meta$`Sample ID` <- NULL

  # Prepare the taxa data
  taxa["LCA"] <- ""
  taxa %>%
    select(ASV, domain, phylum, class, order, family, genus, species, LCA, Contam, ASV_sequence) -> taxa
  taxa           <- as.data.frame(taxa)

  # We need to fill the LCA column starting with 'species', then 'genus', etc
  # until we reach a taxonomy level that isn't 'na' or 'dropped'
  levels <- c("species", "genus", "family", "order", "class", "phylum", "domain")
  for (row in 1:nrow(taxa)) {
    LCA <- NA
    level_ID <- 1
    while(LCA == "dropped" | is.na(LCA)) {
      LCA <- taxa[row, levels[level_ID]]
      level_ID <- level_ID + 1
    }
    taxa[row, "LCA"] <- LCA
  }

  rownames(taxa) <- taxa$ASV
  taxa$ASV       <- NULL
  taxa           <- as.matrix(taxa)

  # Prepare otu data
  otu           <- as.data.frame(otu)
  rownames(otu) <- otu$ASV
  otu[,c('ASV', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'Contam', 'ASV_sequence')] <- list(NULL)

  # Prepare the phylo tree data
  DNA_set = DNAStringSet(seq_tab$ASV_sequence)
  names(DNA_set) = paste0(seq_tab$ASV)

  alignment = AlignSeqs(DNA_set, anchor=NA, processors=cores)

  phang_align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang_align)
  treeNJ <- NJ(dm)

  fit = pml(treeNJ, data=phang_align)
  fitGTR <- update(fit, k=4, inv=0.2)

  if (optimise) {
    # Use 'optim.pml()' to optimise the model parameters
    fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
  }

  # Create the phyloseq object
  OTU = otu_table(otu, taxa_are_rows = TRUE)
  TAX = tax_table(taxa)
  META = sample_data(meta)
  TREE = phy_tree(fitGTR$tree)

  physeq = phyloseq(OTU, TAX, META, TREE)
  saveRDS(physeq, file = paste0('06-report/',voyage,'_',assay,'_phyloseq.rds'))

  #===========================================================================
}

### Done!
print(paste0(voyage, " ", assay, " ", option, " phyloseq object successfully created!"))