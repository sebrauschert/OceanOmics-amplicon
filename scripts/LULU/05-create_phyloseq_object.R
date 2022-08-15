#===========================================================================
# Create phyloseq object (with phyloseq)
# based on here: https://joey711.github.io/phyloseq/import-data.html
# Seb Rauschert
# 22/07/2022
# Modified: Jessica Pearce
# 28/07/2022
#===========================================================================

library(phyloseq)
library(readr)
library(getopt)

# specify input options
# Column 1: the long flag name.  A multi-\link{character} string.
# Column 2: short flag alias of Column 1.  A single-\link{character} string.
# Column 3: 0 = no argument; 1 = required argument; 2 = optional argument
# Column 4: Data type of flag
spec = matrix(c(
  'voyage', 'v', 1, "character",
  'assay' , 'a', 1, "character"), 
  byrow=TRUE, ncol=4)
opt = getopt(spec)

voyage = opt$voyage
assay = opt$assay

# Specify paths to data 
otu_mat    <- paste0('06-lulu/LULU_curated_counts_',voyage,'_',assay,'.csv')
tax_mat    <- paste0('06-lulu/LULU_curated_counts_',voyage,'_',assay,'.csv')
samples_df <- paste0("05-report/",voyage,"_metadata.csv")

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
  select(ASV, kingdom, phylum, class, order, family, genus, LCA, ASV_sequence, Contam) -> taxa
taxa           <- as.data.frame(taxa)
rownames(taxa) <- taxa$ASV
taxa$ASV       <- NULL
taxa           <- as.matrix(taxa)

# Prepare otu data
otu           <- as.data.frame(otu)
rownames(otu) <- otu$ASV
otu[,c('ASV', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'LCA', 'ASV_sequence', 'Contam')] <- list(NULL)

# Create the phyloseq object
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
META = sample_data(meta)

physeq = phyloseq(OTU, TAX, META)
saveRDS(physeq, file = paste0('05-report/',voyage,'_',assay,'_LULU_phyloseq.rds'))

#===========================================================================
