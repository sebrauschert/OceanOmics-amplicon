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