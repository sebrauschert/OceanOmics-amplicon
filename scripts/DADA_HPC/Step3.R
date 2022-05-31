library(dada2)
# Load all Rdata files from step 2 into a list, where the key is the sample name and the value is the output of mergePairs() (`merger`)
path <- '.'
files_to_merge <- list.files(path, pattern='*merged.Rdata')

mergers <- list()
for( f in files_to_merge) {
    sample <- gsub('merged.Rdata', '', f)
    load(f) # now we have a new object called merger
    mergers[[sample]] <- merger
}

## Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

write.csv(seqtab, 'Seqtab.csv', row.names=FALSE)
seqtab <- read.csv('Seqtab.csv', row.names=NULL)
seq_table_nochim <- removeBimeraDenovo(seqtab, 
                                        method = "pooled", 
                                        multithread = 24, 
                                        verbose = TRUE)

asv_seqs <- colnames(seq_table_nochim)
asv_headers <- vector(dim(seq_table_nochim)[2], mode="character")
  
for (i in 1:dim(seq_table_nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, 'asv.fasta') ## input for blastn
  
