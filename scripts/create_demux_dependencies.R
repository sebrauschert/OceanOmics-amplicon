suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# Define options for command line
option_list = list(
  make_option(c("-v", "--voyage"), action="store", default=NA, type='character',
              help="voyage identifier code"),
  make_option(c("-a", "--assay"), action="store", default=NA, type='character',
              help="assay, e.g. '16S' or 'MiFish"))

opt = parse_args(OptionParser(option_list=option_list))

voyage   <- opt$voyage
assay    <- opt$assay


# All we need for now is the sample names and the forward and reverse
barcodes_all <- read_csv(paste0('00-raw-data/indices/', voyage, '_indices.csv')) 

#==========================================================================
# index forward and reverse files are not correct in the original OneDrive location
# here we reformat the files based on the ringtest_indices.csv file
# Seb Rauschert
#==========================================================================

# We loop over two assays and extract the sample Fw and Rv indices and save them in separate fasta files of the format:
# >Fw_name
# Sequence
# >Rv_name
# Sequence
a <- assay
barcodes <- barcodes_all %>%
  filter(assay %in% a)

# Replace special characters with underscores
barcodes$sampleID <- gsub("[^[:alnum:]_]", "_", barcodes$sampleID)

# Append '_1', '_2', etc. to any duplicate samples
duplicate_rows <- duplicated(barcodes$sampleID) | duplicated(barcodes$sampleID, fromLast = TRUE)
underscore_numbers <- ave(seq_along(barcodes$sampleID), barcodes$sampleID, FUN = function(x) seq_along(x))
barcodes$sampleID[duplicate_rows] <- paste0(barcodes$sampleID[duplicate_rows], "_", underscore_numbers[duplicate_rows])

# This will generate a .fa file that searches for both the Fw and the Rv file in R1; whilst keeping the same sample name.
cat(paste(paste0('>',assay_name,'_',barcodes$Fw),
          barcodes$index_seq_Fw,
          paste0('>',assay_name,'_',barcodes$Rv),
          barcodes$index_seq_Rv, sep='\n'), sep = '\n' , file = paste0('00-raw-data/indices/', voyage, '_', assay_name, '_Fw.fa'))

cat(paste(paste0('>',assay_name,'_',barcodes$Rv),
          barcodes$index_seq_Rv,
          paste0('>',assay_name,'_',barcodes$Fw),
          barcodes$index_seq_Fw, sep='\n'), sep = '\n' , file = paste0('00-raw-data/indices/', voyage, '_', assay_name, '_Rv.fa'))

# This section creates a file to rename the demultiplexed files to reflect the sample name, including the assay
cat(paste0(assay_name,'_', barcodes$Fw, '-',assay_name,"_", barcodes$Rv, '.R[12].fq.gz ', barcodes$sampleID ,'_',assay_name,'.#1.fq.gz'),
    paste0(assay_name,'_', barcodes$Rv, '-',assay_name,"_", barcodes$Fw, '.R[12].fq.gz ', barcodes$sampleID ,'_',assay_name,'_reverse.#1.fq.gz'),
    sep='\n',
    file=paste0('00-raw-data/indices/Sample_name_rename_pattern_',voyage,'_',assay_name, '.txt'))
