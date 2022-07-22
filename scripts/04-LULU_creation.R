#-----------------------------------------------------------------------------------------------------------------------------------
# LULU
# https://github.com/tobiasgf/lulu
#-----------------------------------------------------------------------------------------------------------------------------------

library(lulu)
library(readr)

#-----------------------------------------------------------------------------------------------------------------------------------
# 16S
otutab <- read.csv("04-taxa/RS_16S_tax_table_taxonkit.csv",sep=',',header=TRUE,as.is=TRUE, row.names = 1)
otutab[, c('taxid', 'kindom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'LCA', "ASV_sequence",  "Contam")] <- list(NULL)

matchlist_name <- read.table('06-lulu/RS_16S_match_list.txt', header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

curated_result <- lulu(otutab, matchlist_name)
saveRDS(curated_result, file = '06-lulu/16S_LULU_object.rds')

#-----------------------------------------------------------------------------------------------------------------------------------
# MiFish
otutab <- read.csv("04-taxa/RS_MiFish_tax_table_taxonkit.csv",sep=',',header=TRUE,as.is=TRUE, row.names = 1)
otutab[, c('taxid', 'kindom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'LCA', "ASV_sequence",  "Contam")] <- list(NULL)

matchlist_name <- read.table('06-lulu/RS_MiFish_match_list.txt', header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

curated_result <- lulu(otutab, matchlist_name)
saveRDS(curated_result, file = '06-lulu/MiFish_LULU_object.rds')
#-----------------------------------------------------------------------------------------------------------------------------------

rm(list =ls())

LULU_16S    <- readRDS('06-lulu/16S_LULU_object.rds')
LULU_MiFish <- readRDS('06-lulu/MiFish_LULU_object.rds')

RS_16S_tax <- read_csv('04-taxa/RS_16S_tax_table_taxonkit.csv') %>%
  select(ASV, taxid, kindom, phylum, class, order, family, genus, LCA, ASV_sequence, Contam)

RS_MiFish_tax <- read_csv('04-taxa/RS_MiFish_tax_table_taxonkit.csv') %>%
  select(ASV, taxid, kindom, phylum, class, order, family, genus, LCA, ASV_sequence, Contam)


LULU_16S    <- as_tibble(tibble::rownames_to_column(LULU_16S$curated_table, "ASV"))
LULU_MiFish <- as_tibble(tibble::rownames_to_column(LULU_MiFish$curated_table, "ASV"))

LULU_curated_16S    <- as_tibble(merge(RS_16S_tax, LULU_16S, by = 'ASV'))
LULU_curated_MiFish <- as_tibble(merge(RS_MiFish_tax, LULU_MiFish, by = 'ASV'))

write_csv(LULU_curated_16S, '06-lulu/LULU_curated_16S_counts.csv')
write_csv(LULU_curated_MiFish, '06-lulu/LULU_curated_MiFish_counts.csv')

