# Merge the lineage and LCA table

library(tidyverse)
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

Decontam_tab <- read_csv(paste0("04-taxa/",voyage,"_",assay,"_decontam_table.csv"))
Lineage_tab <- read_tsv(paste0("04-taxa/taxonkit_lineage_",voyage,"_",assay,".tsv"))
  
tab_names <- names(Decontam_tab)[4:length(names(Decontam_tab))]
  
as_tibble(cbind(Lineage_tab[,2:length(names(Lineage_tab))], Decontam_tab)) %>%  
  select(c(ASV, names(Lineage_tab), LCA, tab_names)) %>%
  arrange(ASV) %>%
  write_csv(paste0("04-taxa/",voyage,"_",assay,"_tax_table_taxonkit.csv"))
