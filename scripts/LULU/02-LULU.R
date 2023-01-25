#-----------------------------------------------------------------------------------------------------------------------------------
# LULU
# https://github.com/tobiasgf/lulu
#-----------------------------------------------------------------------------------------------------------------------------------

library(lulu)
library(readr)
library(getopt)

# Set working directory for this script
# this is necessary for the docker version of this script
if(Sys.getenv("ANALYSIS") == ""){

  next

}else{

  setwd(Sys.getenv("ANALYSIbS"))

}


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

otutab <- read.csv(paste0("04-LULU/",voyage, "_", assay, "_lulu_table.csv"),sep=',',header=TRUE,as.is=TRUE, row.names = 1, check.names = FALSE)
otutab[, c("ASV_sequence")] <- list(NULL)

matchlist_name <- read.table(paste0('04-LULU/',voyage,'_',assay,'_match_list.txt'), header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

curated_result <- lulu(otutab, matchlist_name)
saveRDS(curated_result, file = paste0('04-LULU/',voyage,'_',assay,'_LULU_object.rds'))

LULU_tab    <- as_tibble(tibble::rownames_to_column(curated_result$curated_table, "ASV"))
write_csv(LULU_tab, paste0('04-LULU/LULU_curated_counts_',voyage,'_',assay,'.csv'))
