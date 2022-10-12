#-----------------------------------------------------------------------------------------------------------------------------------
# LULU
# https://github.com/tobiasgf/lulu
#-----------------------------------------------------------------------------------------------------------------------------------

library(lulu)
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

otutab <- read.csv(paste0("04-taxa/",voyage,"_",assay,"_tax_table_taxonkit.csv"),sep=',',header=TRUE,as.is=TRUE, row.names = 1, check.names = FALSE)
otutab[, c('taxid', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'LCA', "ASV_sequence",  "Contam")] <- list(NULL)
  
matchlist_name <- read.table(paste0('06-lulu/',voyage,'_',assay,'_match_list.txt'), header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
  
curated_result <- lulu(otutab, matchlist_name)
saveRDS(curated_result, file = paste0('06-lulu/',voyage,'_',assay,'_LULU_object.rds'))
  
LULU_tab    <- readRDS(paste0('06-lulu/',voyage,'_',assay,'_LULU_object.rds'))

tax_tab <- read_csv(paste0('04-taxa/',voyage,'_',assay,'_tax_table_taxonkit.csv')) %>%
  select(ASV, taxid, kingdom, phylum, class, order, family, genus, LCA, ASV_sequence, Contam)
  
LULU_tab    <- as_tibble(tibble::rownames_to_column(LULU_tab$curated_table, "ASV"))
LULU_curated    <- as_tibble(merge(tax_tab, LULU_tab, by = 'ASV'))
write_csv(LULU_curated, paste0('06-lulu/LULU_curated_counts_',voyage,'_',assay,'.csv'))

