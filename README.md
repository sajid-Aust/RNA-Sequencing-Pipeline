# RNA-Seq Analysis Pipeline: From Raw Reads to Gene Counts
This repository contains a command-line pipeline for analyzing RNA sequencing (RNA-Seq) data. The pipeline uses HISAT2 for aligning reads to a reference genome and HTSeq-count for quantifying gene expression levels.

## Pipeline Steps:

    1. Create Directories: mkdir -p sam_file bam_file count_file genome_index: Sets up the necessary folders to              organize the output files at different stages of the pipeline.
