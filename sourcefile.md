# To download the data:

# Create environment
conda create -n Viral_triggers_of_Type_1_diabetes -c bioconda -c conda-forge \
  sra-tools fastqc multiqc hisat2 samtools trimmomatic subread -y
conda activate Viral_triggers_of_Type_1_diabetes

# Folder setup
mkdir -p ~/0_Viral_triggers_of_Type_1_diabetes/{data,fastq,trimmed,aligned,counts,logs,qc}
cd ~/0_Viral_triggers_of_Type_1_diabetes/data

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879333/SRR30879333.lite.1 # SRX26276619: EndoC-βH1, Control, biol rep1, unifected

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879324/SRR30879324.lite.1) # SRX26276620: EndoC-βH1, Control, biol rep2, unifected

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879330/SRR30879330.lite.1 # SRX26276622: EndoC-βH1, CVB4-E2 infected, biol rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879321/SRR30879329.lite.1 # SRX26276629: EndoC-βH1, CVB4-E2 infected, biol rep2

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879327/SRR30879327.lite.1 # SRX26276628: EndoC-βH1, CVB4-JVB infected, biol rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879326/SRR30879326.lite.1) # SRX26276626: EndoC-βH1, CVB4-JVB infected, biol rep2

for r in "${SRR[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done
