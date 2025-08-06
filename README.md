Senescence RNA-seq — Preprocessing & Analysis

This repository contains an end-to-end workflow for RNA-seq analysis of early passage and senescent BJ fibroblasts in the absence or continued presence of doxycycline (20nm).
It includes two main components:

Preprocessing (preprocess.sh)
Downloads SRA runs, builds 4 sample FASTQs, performs QC, alignment, and outputs gene transcription level.

Analysis (SenescenceRNAseq_analysis.Rmd)
Performs differential expression analysis using edgeR, generates PCA and volcano plots, and visualizes gene expression levels.

Data

SRA BioProject: 35933466
BJ fibroblasts transduced with doxycycline inducible empty vector control retroviruses were grown in the absence
or continued presence of doxycycline (20nm), and young-quiescent or senescent cells were collected for RNA-seq. 
Senescent and young doxycycline treated samples were compared for changes in genes expression, and senescent and 
young doxycycline untreated samples were compared separately.

Sample 1: GSM4395596 (SRX7865899)
cell line: BJ Fibroblast
Stage: Young - Quiescent
treatment: Doxycycline Minus

Sample 2: GSM4395597 (SRX7865900)
cell line: BJ Fibroblast
Stage: Young - Quiescent
treatment: Doxycycline Plus

Sample 3: GSM4395598 (SRX7865901)
cell line: BJ Fibroblast
Stage: Senescent
treatment: Doxycycline Minus

Sample 4: GSM4395599 (SRX7865902)
cell line: BJ Fibroblast
Stage: Senescent
treatment: Doxycycline Plus

Reference files:

Transcript annotation: gencode.v48.annotation.gtf
Alignment index: hg38
Related publication:
PMID: 40493192

Directory Output

The preprocessing script creates a structured working directory data_pre_processing/ containing:

raw/      # downloaded and processed SRA FASTQs  
fastq/    # concatenated FASTQs (one per sample)  
aligned/  # STAR-aligned BAMs  
counts/   # gene-level count matrix  
logs/     # logs from trimming, alignment, counting  
qc/       # FastQC output  
STAR_index/  # reference index for alignment  
Notes
