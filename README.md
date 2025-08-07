Viral triggers of Type 1 diabetes — Preprocessing & Analysis

This repository contains an end-to-end workflow for small RNA-sequencing to show that infection of immortalized trophoblast cells with coxsackievirus caused differential regulation of several miRNAs.
It includes two main components:

Preprocessing (preprocess.sh)
Downloads SRA runs, builds 4 sample FASTQs, performs QC, alignment, and outputs gene transcription level.

Analysis (SmallRNAseq_analysis.Rmd)
Performs differential expression analysis using edgeR, generates PCA and volcano plots, and visualizes gene expression levels.

Data

SRA BioProject: 39697323
A novel microRNA promotes coxsackievirus B4 infection of pancreatic β cells

Sample 1: GSM8554348 (SRX26276619)
cell line: EndoC-BetaH1
cell type: Pancreatic beta cells
genotype: Wild type
treatment: Control

Sample 2: GSM8554351 (SRX26276622)
cell line: EndoC-BetaH1
cell type: Pancreatic beta cells
genotype: Wild type
treatment: CVB4-E2 infected

Sample 3: GSM8554354 (SRX26276625)
cell line: EndoC-BetaH1
cell type: Pancreatic beta cells
genotype: Wild type
treatment: CVB4-JVB infected

Sample 4: GSM8554357 (SRX26276628)
cell line: Sw.71
cell type: Trophoblast cells
genotype: Wild type
treatment: Control

Sample 5: GSM8554360 (SRX26276631)
cell line: Sw.71
cell type: Trophoblast cells
genotype: Wild type
treatment: CVB4-JVB infected

Reference files:

Transcript annotation: gencode.v48.annotation.gtf
Alignment index: hg38
Related publication:
PMID: 39697323

Directory Output

The preprocessing script creates a structured working directory data_pre_processing/ containing:

raw/      # downloaded and processed SRA FASTQs  
fastq/    # concatenated FASTQs (one per sample)  
aligned/  # STAR-aligned BAMs  
counts/   # gene-level count matrix  
logs/     # logs from trimming, alignment, counting  
qc/       # FastQC output  
STAR_index/  # reference index for alignment  

