#!/bin/bash

# RNA-seq mini-pipeline: QC → alignment → counting

# -------------------- Setup --------------------


# Define project structure relative to current location
PROJECT_DIR="./data_pre_processing"
mkdir -p "${PROJECT_DIR}"/{raw,fastq,aligned,counts,logs,qc,STAR_index}


cd "${PROJECT_DIR}/raw"


# Group SRA run IDs by biological sample 
EndoC-βH1_control=(SRX26276619)   # SRX26276619
CVB4-E2_in=(SRX26276622)   # SRX26276621
EndoC-βH1_CVB4-JVB_in=(SRX26276625)    # SRX26276625
Sw.71_control=(SRX26276628)    # SRX26276628
Sw.71_CVB4-JVB_in=(SRX26276631)   # SRX26276631

# -------------------- Download & Convert --------------------

# Download .sra files
for r in "${EndoC-βH1_control[@]}" "${CVB4-E2_in[@]}" "${EndoC-βH1_CVB4-JVB_in[@]}" "${Sw.71_control[@]}" "${Sw.71_CVB4-JVB_in[@]}"}"; do
  prefetch "$r"
done

# Convert to gzipped FASTQ

for r in "${EndoC-βH1_control[@]}" "${CVB4-E2_in[@]}" "${EndoC-βH1_CVB4-JVB_in[@]}" "${Sw.71_control[@]}" "${Sw.71_CVB4-JVB_in[@]}"}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done

# Concatenate per-sample FASTQs
cat "${EndoC-βH1_control[@]/%/.fastq.gz}"  > EndoC-βH1_control.fastq.gz
cat "${CVB4-E2_in[@]/%/.fastq.gz}"  > CVB4-E2_in.fastq.gz
cat "${EndoC-βH1_CVB4-JVB_in[@]/%/.fastq.gz}" > EndoC-βH1_CVB4-JVB_in.fastq.gz
cat "${Sw.71_control[@]/%/.fastq.gz}" > Sw.71_control.fastq.gz
cat "${Sw.71_CVB4-JVB_in[@]/%/.fastq.gz}"  > Sw.71_CVB4-JVB_in.fastq.gz

# Move to fastq/ folder
mv CVB4.fastq.gz Endo*.fastq.gz Sw*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ../fastq
fastqc EndoC-βH1_control.fastq.gz CVB4-E2_in.fastq.gz EndoC-βH1_CVB4-JVB_in.fastq.gz Sw.71_control.fastq.gz Sw.71_CVB4-JVB_in.fastq.gz \
  -o ../qc --threads 16

# -------------------- Alignment (STAR) --------------------

cd ../STAR_index
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
unzip GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v46.annotation.gtf.gz
unzip gencode.v46.annotation.gtf.gz

module load star
GENOMEDIR="./RNAseq/genome/"
mdkir -p $GENOMEDIR/STAR
STAR --runThreadN 23 --runMode genomeGenerate --genomeDir $GENOMEDIR/STAR --genomeFastaFiles $GENOMEDIR/GRCh38.primary_assembly.genome.fa --sjdbGTFfile $GENOMEDIR/gencode.v29.primary_assembly.annotation.gtf
cd ../trimmed
STAR --genomeDir indexes/chr10 \
      --readFilesIn EndoC-βH1_control.fastq.gz CVB4-E2_in.fastq.gz EndoC-βH1_CVB4-JVB_in.fastq.gz Sw.71_control.fastq.gz Sw.71_CVB4-JVB_in.fastq.gz \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix alignments/

# -------------------- Quantification (featureCounts) --------------------

cd ..
curl -L -o gencode.v46.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v46.annotation.gtf.gz
gunzip -f gencode.v46.annotation.gtf.gz

featureCounts -T 16 -t exon -g gene_name \
  -a gencode.v46.annotation.gtf \
  -o counts/raw_counts_gene_sym.txt aligned/*.bam \
  &> logs/featureCounts_gene_sym.log

# Format counts matrix
{ printf "GeneSymbol\t"; head -n 2 counts/raw_counts_gene_sym.txt | tail -n 1 | cut -f7-; } > counts/final_counts_symbols.tsv
tail -n +3 counts/raw_counts_gene_sym.txt | \
  awk -v OFS="\t" '{ out=$1; for(i=7;i<=NF;i++) out=out OFS $i; print out }' >> Processed_counts/final_counts_symbols.tsv

sed -i '' '1 s|aligned/||g; 1 s|\.bam||g' counts/final_counts_symbols.tsv

# Done
echo "Pipeline complete. Output saved in: ${PROJECT_DIR}/counts/final_counts_symbols.tsv"
