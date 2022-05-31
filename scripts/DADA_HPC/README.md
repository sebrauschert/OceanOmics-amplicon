# DADA2 pipeline for HPC

With many samples, the standard DADA workflow takes a very long time (days to weeks). Yet some steps of the DADA2 pipeline are highly parellelisable. This workflow is a replacement of the 04-DADA2.R script in cases where there are many samples (>100).
This workflow is based on https://benjjneb.github.io/dada2/bigdata.html but splits up the per-sample processing into separate scripts for ease of submission to SLURM.

In the 'default' DADA2 pipeline, each sample is processed sequentially, then all samples are merged at the end. The workflow here processes the samples simultaneously on a cluster. Step 1 and Step 3 are run on one node, Step 2 is run on many nodes.

What it does:

* Step1.R : learn the error rate from a few samples and 1e8 bases (default DADA2 settings). Write MANY R-scripts, one for each sample. Run this once.
* Step2.SAMPLENAME.R : Each R-script can be submitted to SLURM separately, one script per sample, which runs the actual dada and mergeSamples commands. Run this once for each sample.
* Step3.R : this script loads all the dada2 results from each sample and calculates the final ASVs, writes them to a table and then a fasta file.


# Running the pipeline

Make sure to adjust the `path` variable in all scripts.

Run step 1, wait a few minutes, then you will have many R-scripts in the current directory. Use a SLURM array job or Snakemake/Nextflow to submit them all to the cluster, wait. It should take less than 30 minutes per sample. Then run step 3 and receive the final fasta.


# TODO

* make the path and the file endings more dynamic
* parallelise the filterAndTrim() step in Step1 the same way, which would mean another set of R-scripts
* make all the scripts in Step 2 nicer. Perhaps in their own folder?
