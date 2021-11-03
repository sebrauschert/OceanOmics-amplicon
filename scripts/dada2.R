#'''''''''''''''''''''''''''''''''''''''''''''''''''
#'SETUP FOR DADA2 ANALYSIS
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Install DADA2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.13")

# Install tidyverse
# install.packages('tidyverse')


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'ANALYSIS
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

library(dada2)
library(tidyverse)
library(RColorBrewer)

# Define fish or elasmo

species="elasmo"


path         <- paste0("02-demultiplexed/", species)
raw_forward  <- sort(list.files(path,
                                pattern = ".R1.fq.gz",
                                full.names = TRUE))
sample_names <- sapply(strsplit(basename(raw_forward), "-"), `[`, 1)

# visualise the quality of the reads
plotQualityProfile(raw_forward[1:12])

# place filtered files in filtered/ subdirectory
filtered_path    <- file.path("02-filtered")
filtered_forward <- file.path(filtered_path,
                              paste0(sample_names, "-R1_trimmed.fq.gz"))

# quality filtering
# remove barcode + primers from beginning and end of the reads
# use number of bases to remove from left and from right

# check length of barcodes
tags <- read_csv(paste0("00-raw-data/adapters/sequence_tags_",species,".csv"))
head(tags)
nchar(tags$Sequence) # all 8 bp
len_barcode <- unique(nchar(tags$Sequence))

# check length of primers
primers <- read.csv("00-raw-data/adapters/sequence_primers_fish_elasmo.csv")
head(primers)

# Special case as both fish and elasmo are first letter upper case in the table
selection = stringr::str_to_title(species)

primers <- primers %>%
  filter(Target == selection) %>%
  mutate(len_primer = nchar(Sequence))

len_Fw <- primers$len_primer[1]
len_Rv <- primers$len_primer[2]

# combine length barcode + primer
trim_len_Fw <- len_barcode + len_Fw
trim_len_Rv <- len_barcode + len_Rv

# Trimming based on quality + barcode/primer removal
out <- filterAndTrim(fwd = raw_forward, 
                     filt = filtered_forward, 
                     maxN = 0, 
                     trimLeft = trim_len_Fw, 
                     trimRight = trim_len_Rv, 
                     maxEE = 2, 
                     rm.phix = TRUE, 
                     compress = TRUE, 
                     multithread = TRUE, 
                     verbose = TRUE)

head(out)
# reads.in reads.out
# 16S2P4B-fish.R1.fq.gz      110        89
# BY1_R1-fish.R1.fq.gz    136501    132777
# BY1_R2-fish.R1.fq.gz    103384    100127
# BY1_R3-fish.R1.fq.gz     97543     94420
# BY1_R4-fish.R1.fq.gz    100832     97383
# BY1_R5-fish.R1.fq.gz    113757    110345

# visualise the quality of the trimmed reads
plotQualityProfile(filtered_forward[1:12])

# run seqkit stats on terminal 
# seqkit stats -T /media/Priscila_ext_HD/Data/amplicon/MiSeq_voyage1/fish/02_Filtered/*_trimmed.fq.gz > /media/Priscila_ext_HD/Data/amplicon/MiSeq_voyage1/fish/02_Filtered/trimmed_seqkit_stats.txt 

# Inspect Stats results and compile number into Excel spreadsheet on ShareFile - Minderoo_tags_Voyage1_sequencing

# Learn the error rates
errors_forward <- learnErrors(filtered_forward, multithread = TRUE)

# visualise the estimated error rates
plotErrors(errors_forward, nominalQ = TRUE)

# De-replication
derep_forward        <- derepFastq(filtered_forward, verbose = TRUE)
names(derep_forward) <- sample_names
head(derep_forward)

# Sample inference¶
dada_forward <- dada(derep_forward, 
                     err = errors_forward, 
                     pool = TRUE, 
                     multithread = TRUE, 
                     verbose = TRUE)

# inspect the dada-class object
dada_forward[[1]]

# Construct Sequence Table
seq_table <- makeSequenceTable(dada_forward)
dim(seq_table) #  92 4290

# inspect distribution of sequence lengths
table(nchar(getSequences(seq_table)))

# Remove Chimeras
# if pooling for denoising, should also pool for chimera removal
seq_table_nochim <- removeBimeraDenovo(seq_table, 
                                       method = "pooled", 
                                       multithread = TRUE, 
                                       verbose = TRUE)

# Identified 3133 bimeras out of 4290 input sequences.

dim(seq_table_nochim) # 92 1157

# which percentage of our reads did we keep?
sum(seq_table_nochim) / sum(seq_table)
# 97% of reads kept

dim(seq_table_nochim) [2] / dim(seq_table)[2]
# 27% of ASVs kept

# Overview of counts throughout
get_n <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dada_forward, get_n), rowSums(seq_table), rowSums(seq_table_nochim))  %>%
  as.data.frame() %>%
  mutate(final_perc_reads_retained = round(rowSums(seq_table_nochim)/out[,1]*100, 1))

colnames(track) <- c('input', 'filtered', 'denoised', 'tabled', 'nonchim', 'final_perc_reads_retained')

rownames(track) <- sample_names
head(track)
tail(track)

summary(track$nonchim)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     25   69123   85764   82816  104328  134240 

# plot out tracking of sample reads through stages ####
samps <- row.names(track)
track_df <- data.frame(track) %>%
  mutate(samps = samps,
         site = sapply(strsplit(samps, "_"), `[`, 1)) %>%
  gather('stage', 'reads', c(input, filtered, denoised, tabled, nonchim))

track_boxplot <- ggplot(track_df, aes(forcats::fct_relevel(stage, c('input', 'filtered', 'denoised', 'tabled', 'nonchim')), reads)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = site), position = position_jitter(width = 0.2, height = 0), shape = 21, alpha = 0.7, size = 2) +
  scale_fill_manual("Site", values = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(track_df$site)))) +
  ylab('Number of reads') +
  xlab('Sequencing stage') +
  theme_bw() +
  theme(legend.position = "bottom")


#ggsave(plot = track_boxplot, filename = "04-dada2/figures/track_reads_fish.pdf", height = 5, width = 7)

#Check if all sequences are the same length
plyr::count(unlist(lapply(colnames(seq_table_nochim), function(x) stringi::stri_length(x))))



# Save the ASV sequences as .fa file

asv_seqs <- colnames(seq_table_nochim)
asv_headers <- vector(dim(seq_table_nochim)[2], mode="character")

for (i in 1:dim(seq_table_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0("/data/voyage1",species,".fa"))

# BLAST it 
#blastn -db databases/nt \
#-query voyage1.fa \
#-evalue 1e-3 \
#-num_threads 16 -html > voyage1_results
