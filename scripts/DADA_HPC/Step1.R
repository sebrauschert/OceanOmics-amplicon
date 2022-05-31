library(dada2); packageVersion("dada2")
# CHANGE THE PATH
path <- '.'
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
print(sample.names)

png('before_qualR1.png')
plotQualityProfile(fnFs[1:2])
dev.off()


png('before_qualR2.png')
plotQualityProfile(fnRs[1:2])
dev.off()


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(224,236),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
save(errF, errR, file = 'errFR.Rdata')

png('estimated_errors.png')
plotErrors(errF, nominalQ=TRUE)
dev.off()

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
# Now write MANY R-scripts
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  sink(paste0('Step2.', sam, '.clean.R'))
  cat("library(dada2)\n")
  cat("path <- '/scratch/pawsey0149/pbayer/amplicon_fun/OceanOmics-amplicon/taraAmplicon_Amplicon_pbayer/PRJEB6610/results/fastq/'\n")
  cat("fnFs <- sort(list.files(path, pattern='_1.fastq.gz', full.names = TRUE))\n")
  cat("fnRs <- sort(list.files(path, pattern='_2.fastq.gz', full.names = TRUE))\n")
  cat("sample.names <- sapply(strsplit(basename(fnFs), '_'), `[`, 1)\n")
  cat("filtFs <- file.path(path, 'filtered', paste0(sample.names, '_F_filt.fastq.gz'))\n")
  cat("filtRs <- file.path(path, 'filtered', paste0(sample.names, '_R_filt.fastq.gz'))\n")
  cat("names(filtFs) <- sample.names\n")
  cat("names(filtRs) <- sample.names\n")
  cat("load('errFR.Rdata')\n")
  cat("derepF <- derepFastq(filtFs[['", sam, "']])\n", sep='')
  cat("ddF <- dada(derepF, err=errF, multithread=TRUE)\n")
  cat("derepR <- derepFastq(filtRs[['", sam, "']])\n", sep='')
  cat("ddR <- dada(derepR, err=errR, multithread=TRUE)\n")
  cat("merger <- mergePairs(ddF, derepF, ddR, derepR)\n")
  cat("save(merger, file='", paste0(sam, "merged.Rdata"), "')\n", sep='')
  sink()
}

# now run all these scripts via SLURM or Snakemake or etc.
