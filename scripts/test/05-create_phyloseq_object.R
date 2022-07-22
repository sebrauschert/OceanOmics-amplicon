#===========================================================================
# Create phyloseq object (with phyloseq)
# based on here: https://joey711.github.io/phyloseq/import-data.html
# Seb Rauschert
# 22/07/2022
#===========================================================================

# Specify paths to data 
otu_mat    <- '06-lulu/LULU_curated_16S_counts.csv'
tax_mat    <- '06-lulu/LULU_curated_16S_counts.csv'
samples_df <- "05-report/RS_metadata.csv"

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
  select(ASV, taxid, kindom, phylum, class, order, family, genus, LCA, ASV_sequence, Contam) -> taxa
taxa           <- as.data.frame(taxa)
rownames(taxa) <- taxa$ASV
taxa$ASV       <- NULL
taxa           <- as.matrix(taxa)

# Prepare otu data
otu           <- as.data.frame(otu)
rownames(otu) <- otu$ASV
otu[,c('ASV', 'taxid', 'kindom', 'phylum', 'class', 'order', 'family', 'genus', 'LCA', 'ASV_sequence', 'Contam')] <- list(NULL)

# Create the phyloseq object
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
META = sample_data(meta)

physeq = phyloseq(OTU, TAX, META)
saveRDS(physeq, file = '05-report/RS_16S_LULU_phyloseq.rds')

#===========================================================================
# MiFish

# Specify paths to data 
otu_mat    <- '06-lulu/LULU_curated_MiFish_counts.csv'
tax_mat    <- '06-lulu/LULU_curated_MiFish_counts.csv'
samples_df <- "05-report/RS_metadata.csv"

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
  select(ASV, taxid, kindom, phylum, class, order, family, genus, LCA, ASV_sequence, Contam) -> taxa
taxa           <- as.data.frame(taxa)
rownames(taxa) <- taxa$ASV
taxa$ASV       <- NULL
taxa           <- as.matrix(taxa)

# Prepare otu data
otu           <- as.data.frame(otu)
rownames(otu) <- otu$ASV
otu[,c('ASV', 'taxid', 'kindom', 'phylum', 'class', 'order', 'family', 'genus', 'LCA', 'ASV_sequence', 'Contam')] <- list(NULL)

# Create the phyloseq object
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
META = sample_data(meta)

physeq = phyloseq(OTU, TAX, META)
saveRDS(physeq, file = '05-report/RS_MiFish_LULU_phyloseq.rds')
