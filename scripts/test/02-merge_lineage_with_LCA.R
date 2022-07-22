# Merge the lineage and LCA table

library(tidyverse)

RS16S <- read_csv('04-taxa/RS_16S_decontam_table.csv')
RS16S_lineage <- read_tsv('04-taxa/taxonkit_taxa_lineage_16S.tsv')

rs16s_names <- names(RS16S)[4:length(names(RS16S))] 

#as_tibble(merge(RS16S_lineage, RS16S, by ='taxid')) %>%
  
as_tibble(cbind(RS16S_lineage[,2:length(names(RS16S_lineage))], RS16S)) %>%  
  select(c(ASV, names(RS16S_lineage), LCA, rs16s_names)) %>%
  #mutate(ASV = paste0('ASV_', 1:dim(RS16S)[1])) %>%
  arrange(ASV) %>%
  write_csv('04-taxa/RS_16S_tax_table_taxonkit.csv')


# MiFish
RSMiFish <- read_csv('04-taxa/RS_MiFish_decontam_table.csv')
RSMiFish_lineage <- read_tsv('04-taxa/taxonkit_taxa_lineage_MiFish.tsv')

rsMiFish_names <- names(RSMiFish)[4:length(names(RSMiFish))] 

#as_tibble(merge(RSMiFish_lineage, RSMiFish, by ='taxid')) %>%

as_tibble(cbind(RSMiFish_lineage[,2:length(names(RSMiFish_lineage))], RSMiFish)) %>%  
  select(c(ASV, names(RSMiFish_lineage), LCA, rsMiFish_names)) %>%
  #mutate(ASV = paste0('ASV_', 1:dim(RSMiFish)[1])) %>%
  arrange(ASV) %>%
  write_csv('04-taxa/RS_MiFish_tax_table_taxonkit.csv')

rm(list=ls())

RS_16S    <- read_csv('04-taxa/RS_16S_tax_table_taxonkit.csv')
RS_MiFish <- read_csv('04-taxa/RS_MiFish_tax_table_taxonkit.csv')

RS_16S %>%
  select(!c('taxid', 'kindom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'LCA', "ASV_sequence",  "Contam")) %>%
  

as_tibble(cbind(RSMiFish_lineage[,2:length(names(RSMiFish_lineage))], RSMiFish)) %>%  
  select(c(ASV, names(RSMiFish_lineage), LCA, rsMiFish_names)) %>%
  arrange(ASV) -> TEST

