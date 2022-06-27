#........................................................................
# ANALYSIS OF AMPLICON DATA
#
# Data Info: 
# Analyst  : 
# Date     : 
#........................................................................

library(dada2) 
library(tidyverse) 
library(RColorBrewer) 
library(readr)

# DADA2 pipeline as function to 
# enable writing a loop 
dada2_analysis <- function(voyageID = voyageID, 
                           assay = assay, 
                           site = site){
  
  
  # define path
  path         <- paste0(getwd(),"/02-demultiplexed/", voyageID, "/", assay, "/", site)
  list.files(path)
  
  # loading index file
  tags <- read_csv(paste0("00-raw-data/indices/",voyageID,"_indices.csv"))
  
  # read in fastq files
  fnFs <- sort(list.files(path, pattern="1.fq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="2.fq", full.names = TRUE))
  
  # extract the short sample name from the filename
  sample.names_Fs <- as.character(sapply(basename(fnFs), function(x) unlist(stringr::str_remove(x,paste0("_",assay,".1.fq.gz")))))
  sample.names_Rs <- as.character(sapply(basename(fnRs), function(x) unlist(stringr::str_remove(x,paste0("_",assay,".2.fq.gz")))))
  
  # Take a random subset of the samples and save one quality plot at 
  # a time so we can later on add them to the analysis report
  set.seed(4)
  
  for(i in sample(1:length(fnFs), 3, replace=FALSE)){
    
    qualityprofile_Fs <- plotQualityProfile(fnFs[i])
    
    ggsave(plot = qualityprofile_Fs, 
           filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "_qualityprofile_Fs_", i, "_", assay,"_", site,"_raw.png"), 
           height = 5, 
           width = 7)
    
    qualityprofile_Rs <- plotQualityProfile(fnRs[i])
    
    ggsave(plot = qualityprofile_Rs, 
           filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "_qualityprofile_Rs_", i, "_", assay,"_", site, "_raw.png"), height = 5, width = 7)
    
  }
  
  # check barcode lengths
  head(tags)
  if(assay=="MiFish"){
    tags <- subset(tags, assay=="MiFish")
  } else {
    tags <- subset(tags, assay=="16S")
  }
  
  nchar(tags$index_seq_Fw) # all 8 bp
  len_barcode_Fw <- unique(nchar(tags$index_seq_Fw))
  nchar(tags$index_seq_Rv) # all 8 bp
  len_barcode_Rv <- unique(nchar(tags$index_seq_Rv))
  
  nchar(tags$full_primer_seq_Fw)
  len_primer_Fw <- unique(nchar(tags$full_primer_seq_Fw))
  nchar(tags$full_primer_seq_Rv)
  len_primer_Rv <- unique(nchar(tags$full_primer_seq_Rv))
  
  # trim length of primer without barcode
  trim_len_Fw <- len_primer_Fw - len_barcode_Fw
  trim_len_Rv <- len_primer_Rv - len_barcode_Rv
  trim_len_Fw
  trim_len_Rv
  
  # Assigns file names and place filtered files in filtered/sub directory
  filtered_path    <- file.path(paste0(getwd(),"/03-dada2/", voyageID, "/filtered_", voyageID, "_", assay, "/", site))
  filtFs <- file.path(filtered_path,
                      paste0(sample.names_Fs, "_", assay, "_1_trimmed.fq.gz"))
  filtRs <- file.path(filtered_path,
                      paste0(sample.names_Rs, "_", assay, "_2_trimmed.fq.gz"))
  names(filtFs) <- sample.names_Fs
  names(filtRs) <- sample.names_Rs
  
  # Trimming based on quality + barcode/primer removal
  ## the maxEE = expected error - filtering based on quality
  out <- filterAndTrim(fnFs, filtFs, 
                       fnRs, filtRs, 
                       trimLeft = trim_len_Fw,
                       trimRight = trim_len_Rv,
                       minLen = 80,
                       maxN=0, 
                       maxEE=c(2,2), 
                       truncQ=2, 
                       rm.phix=TRUE,
                       compress=TRUE, 
                       multithread=TRUE) # On Windows set multithread=FALSE
  
  head(out)
  #......................................................................................
  # CHECKPOINT Save the result
  saveRDS(out, paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay,"_", site, "_filterAndTrim_out.rds"))
  #out <- readRDS(paste0('03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay,"_", site, "_filterAndTrim_out.rds'))
  #......................................................................................
  
  
  # visualize the quality of the trimmed reads
  # Take a random subset of the samples and save one quality plot at 
  # a time so we can later on add them to the analysis report
  
  set.seed(4)
  
  for(i in sample(1:length(fnFs), 3, replace=FALSE)){
    
    qualityprofile_Fs <- plotQualityProfile(filtFs[i])
                                            
    ggsave(plot = qualityprofile_Fs, 
           filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "_qualityprofile_Fs_", i, "_", assay,"_", site,"_trimmed.png"), 
           height = 5, 
           width = 7)
    
    qualityprofile_Rs <- plotQualityProfile(filtRs[i])
    
    ggsave(plot = qualityprofile_Rs, 
           filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "_qualityprofile_Rs_", i, "_", assay,"_", site, "_trimmed.png"), height = 5, width = 7)
    
  }
  
  
  
  # Learn the error rates
  errors_forward <- learnErrors(filtFs, multithread = TRUE)
  errors_reverse <- learnErrors(filtRs, multithread = TRUE)
  
  #......................................................................................
  # CHECKPOINT Save the result
  save(errors_forward, errors_reverse, file = paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay,"_", site, "_error_rates.RData"))
  #load(paste0("03-dada2", voyageID, "/tmpfiles/", voyageID, "_", assay,"_", site, "_error_rates.RData"))
  #......................................................................................
  
  # visualise the estimated error rates
  errorsplot_Fs <- plotErrors(errors_forward, nominalQ = TRUE)
  ggsave(plot = errorsplot_Fs, 
         filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "_errorsplot_Fs_", assay,"_", site,".png"), 
         height = 5, 
         width = 7)
  errorsplot_Rs <- plotErrors(errors_reverse, nominalQ = TRUE)
  ggsave(plot = qualityprofile_Rs, 
         filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "_errorsplot_Rs_", assay,"_", site, ".png"), 
         height = 5, 
         width = 7)
  
  # De-replication - check if we should do this if we want abundance estimates
  ## supposedly for high resolution sample interference from amplicon data - removes sequencing errors
  derep_forward        <- derepFastq(filtFs, verbose = TRUE)
  names(derep_forward) <- sample.names_Fs
  head(derep_forward)
  
  derep_reverse        <- derepFastq(filtRs, verbose = TRUE)
  names(derep_reverse) <- sample.names_Rs
  head(derep_reverse)
  
  #......................................................................................
  # CHECKPOINT Save the result
  save(derep_forward, derep_reverse, file = paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay,"_", site, "_dereplicated.RData"))
  #load(paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay,"_", site, "_dereplicated.RData"))
  #......................................................................................
  
  # Sample inference
  dada_forward <- dada(derep_forward, 
                       err = errors_forward, 
                       pool = TRUE, 
                       multithread = TRUE, 
                       verbose = TRUE)
  
  dada_reverse <- dada(derep_reverse, 
                       err = errors_reverse, 
                       pool = TRUE, 
                       multithread = TRUE, 
                       verbose = TRUE)
  
  #......................................................................................
  # CHECKPOINT Save the result
  save(dada_forward, dada_reverse, file = paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay, "_", site, "_core_sample_inference.RData"))
  #load(paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay, "_", site, "_core_sample_inference.RData"))
  #......................................................................................
  
  # inspect the dada-class object
  dada_forward[[1]]
  dada_reverse[[1]]
  
  # merge paired end reads
  mergers <- mergePairs(dada_forward, 
                        filtFs, 
                        dada_reverse, 
                        filtRs,
                        minOverlap = 5,
                        verbose=TRUE)
  
  #......................................................................................
  # CHECKPOINT Save the result
  save(mergers, file = paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay, "_", site,"_merged.RData"))
  # load(paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay, "_", site,"_merged.RData"))
  #......................................................................................
  
  
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])
  
  # Construct Sequence Table
  seq_table <- makeSequenceTable(mergers)
  dim(seq_table)
  
  # inspect distribution of sequence lengths
  table(nchar(getSequences(seq_table)))
  
  mean(nchar(getSequences(seq_table)))
  median(nchar(getSequences(seq_table)))
  
  seq_dist <- as.data.frame(nchar(getSequences(seq_table)))
  seq_dist
  
  # Create histogram of sequence length distributions
  seq_hist <- ggplot(seq_dist, aes(nchar(getSequences(seq_table)))) +
    geom_histogram(bins = 100, ) +
    ylab('Number of reads') +
    xlab('Sequence length (bp)') +
    theme(text = element_text(size=25))
  seq_hist
  
  # Save plot
  ggsave(plot = seq_hist,
         filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "_seq_distribution_",assay,"_",site,"_ASVs.png"),
         height = 10,
         width = 12)

  # filter for amplicon length: The 16S and MiFish primers each have a specific range of base pairs
  ### 16S = 178 - 228
  ### MiFish = 163 - 185
  
   if(assay=="MiFish"){
    seq_table2 <- seq_table[,nchar(colnames(seq_table)) %in% 163:185]
  } else {
    seq_table2 <- seq_table[,nchar(colnames(seq_table)) %in% 178:228]
  }
  
  # Remove Chimeras
  # if pooling for denoising, should also pool for chimera removal
  seq_table_nochim <- removeBimeraDenovo(seq_table2, 
                                         method = "pooled", 
                                         multithread = TRUE, 
                                         verbose = TRUE)
  
  dim(seq_table_nochim)
  
  #......................................................................................
  # CHECKPOINT Save the result
  save(seq_table_nochim, file = paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay, "_", site,"_seq_table_nochim.RData"))
  # load(paste0("03-dada2/", voyageID, "/tmpfiles/", voyageID, "_", assay, "_", site,"_seq_table_nochim.RData"))
  #......................................................................................
  
  # which percentage of our reads did we keep?
  sum(seq_table_nochim) / sum(seq_table2)
  dim(seq_table_nochim) [2] / dim(seq_table2)[2]
  
  ## Overview of counts throughout
  get_n <- function(x) sum(getUniques(x))
  
  #forward reads track
  track_Fs <- cbind(out, sapply(dada_forward, get_n), rowSums(seq_table2), rowSums(seq_table_nochim))  %>%
    as.data.frame() %>%
    mutate(final_perc_reads_retained = round(rowSums(seq_table_nochim)/out[,1]*100, 1))
  
  colnames(track_Fs) <- c('input', 'filtered', 'denoised', 'tabled', 'nonchim', 'final_perc_reads_retained')
  rownames(track_Fs) <- sample.names_Fs
  head(track_Fs)
  tail(track_Fs)
  write.table(track_Fs, file = paste0("03-dada2/", voyageID, "/QC_plots/Track_reads_Fw_",assay,"_", site))
  
  #reverse reads track
  track_Rs <- cbind(out, sapply(dada_reverse, get_n), rowSums(seq_table2), rowSums(seq_table_nochim))  %>%
    as.data.frame() %>%
    mutate(final_perc_reads_retained = round(rowSums(seq_table_nochim)/out[,1]*100, 1))
  
  colnames(track_Rs) <- c('input', 'filtered', 'denoised', 'tabled', 'nonchim', 'final_perc_reads_retained')
  rownames(track_Rs) <- sample.names_Rs
  head(track_Rs)
  tail(track_Rs)
  write.table(track_Rs, file = paste0("03-dada2/", voyageID, "/QC_plots/Track_reads_Rs_",assay,"_", site))
  
  summary(track_Fs$nonchim)
  summary(track_Rs$nonchim)
  
  # plot out tracking of sample reads through stages ####
  samps_Fs <- row.names(track_Fs)
  track_df_Fs <- data.frame(track_Fs) %>%
    mutate(samps = samps_Fs,
           site = sapply(strsplit(samps, "_"), `[`, 1)) %>%
    gather('stage', 'reads', c(input, filtered, denoised, tabled, nonchim))
  
  samps_Rs <- row.names(track_Rs)
  track_df_Rs <- data.frame(track_Rs) %>%
    mutate(samps = samps_Rs,
           site = sapply(strsplit(samps, "_"), `[`, 1)) %>%
    gather('stage', 'reads', c(input, filtered, denoised, tabled, nonchim))
  
  # create plots
  track_boxplot_Fw <- ggplot(track_df_Fs, aes(forcats::fct_relevel(stage, c('input', 'filtered', 'denoised', 'tabled', 'nonchim')), reads)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = site), position = position_jitter(width = 0.2, height = 0), shape = 21, alpha = 0.7, size = 2) +
    scale_fill_manual("Site", values = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(track_df_Fs$site)))) +
    ylab('Number of reads') +
    xlab('Sequencing stage') +
    theme_bw() +
    theme(legend.position = "bottom")
  
  track_boxplot_Rv <- ggplot(track_df_Rs, aes(forcats::fct_relevel(stage, c('input', 'filtered', 'denoised', 'tabled', 'nonchim')), reads)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = site), position = position_jitter(width = 0.2, height = 0), shape = 21, alpha = 0.7, size = 2) +
    scale_fill_manual("Site", values = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(track_df_Rs$site)))) +
    ylab('Number of reads') +
    xlab('Sequencing stage') +
    theme_bw() +
    theme(legend.position = "bottom")
  
  
  # save plots
  ggsave(plot = track_boxplot_Fw, 
         filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "Samples_through_stages_Fw_",assay,"_", site,".png"), 
         height = 10, 
         width = 12)
  ggsave(plot = track_boxplot_Rv, 
         filename = paste0("03-dada2/", voyageID, "/QC_plots/", voyageID, "Samples_through_stages_Rv_",assay,"_", site,".png"), 
         height = 10, 
         width = 12)
  
  
  #Check if all sequences are the same length
  plyr::count(unlist(lapply(colnames(seq_table_nochim), function(x) stringi::stri_length(x))))
  
  # Save the ASV sequences as .fa file
  asv_seqs <- colnames(seq_table_nochim)
  asv_headers <- vector(dim(seq_table_nochim)[2], mode="character")
  
  for (i in 1:dim(seq_table_nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
  }
  
  asv_final_table <- seq_table_nochim
  colnames(asv_final_table) <- asv_headers
  IDs <- rownames(asv_final_table)
  
  as_tibble(asv_final_table) %>%
    mutate(sample_id = IDs) %>%
    select(sample_id, asv_headers) -> asv_final_table
  
  
  #--------------------------------------------------------------------------------------------------------------------------------------------------
  # SAVE RESULTS
  
  # making and writing out a fasta of our final ASV seqs:
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  write(asv_fasta, paste0("03-dada2/", voyageID, "/", voyageID, "_",assay,"_", site,".fa")) ## input for blastn
  
  # Prepare ASV table for LCA (Lowest Common Ancestor?)
  ## Here we rename the sequences to ASV with an ID, to match the blast results
  colnames(seq_table_nochim) <- asv_headers
  
  # We need to transpose, so the rows are the sequences, whereas the columns are the IDs
  asv_for_lca <- as.data.frame(t(seq_table_nochim))
  
  # Making sure that we follow the nomenclature for LCA
  headers_lca <- c('ASV', names(asv_for_lca))
  
  # Capture the IDs / sample names
  ID <- rownames(asv_for_lca)
  
  # Execute it all and create the ASV table
  asv_for_lca <- asv_for_lca %>%
    mutate(`ASV`= ID) %>%
    select((headers_lca)) %>%
    as_tibble()
  
  asv_for_lca[,1] <- str_remove(as.vector(unlist(asv_for_lca[,1])) , ">")
  asv_for_lca$ASV_sequence <- asv_seqs
  write_delim(asv_for_lca, paste0("03-dada2/", voyageID, "/",voyageID, "_asv_final_table_",assay,"_", site,".tsv"), delim = '\t')
  
}

#........................................................................
# Running DADA2

voyages = "ABV4"
assays  = c("16S", "MiFish")
sites = c("AUV01_N_10_Wallabi", "AUV01_Pelsaert", "AUV03_N_40_Wallabi", "AUV03_Pelsaert", "AUV04_N_15_Wallabi", "AUV04_Pelsaert", "AUV06_N_35_Wallabi", 
          "AUV06_Pelsaert", "C-153_N_40_Easter", "Controls", "D1_Jurien", "D2_Jurien", "D3_Jurien", "D2_GC", "D5_GC", "Deep3_Wallabi", "DeepChannel_Wallabi", 
          "DeepKelp_Wallabi", "EastWallabi_Wallabi", "E_D1_N_40_Easter", "E_D2_N_40_Easter", "E_S1_N_10_Easter", "E_S2_N_10_Easter", "E_S3_N_10_Easter", 
          "ES_Easter", "Kelp1_Wallabi", "PW_D01_50_Pelsaert", "P_W_D1_N_40_Pelsaert", "PW_D_40_Pelsaert", "P_W_S_10_Pelsaert", "PW_S_10_Pelsaert", 
          "S1_GC", "S2_GC", "S3_GC", "S4_GC", "TurtleBay_Wallabi", "WA115_10_Wallabi", "WA115_Pelsaert", "WA176_Pelsaert", "WA177_10_Wallabi", 
          "WA177_Pelsaert", "WA181_Pelsaert", "WA46_Pelsaert", "WestWallabi_Wallabi")

# After running the function below, this loop will run the full analysis across
# all voyages and assays

for(voyage in voyages) {
  
  for(assay in assays){
    
    for(site in sites){
    
    dada2_analysis(voyage, assay, site)
      
    }
  }
}
