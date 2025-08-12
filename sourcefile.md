To download the data:

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879333/SRR30879333.lite.1 # SRX26276619: EndoC-βH1, Control, biol rep1, unifected

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879330/SRR30879330.lite.1 # SRX26276622: EndoC-βH1, CVB4-E2 infected, biol rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879327/SRR30879327.lite.1 # SRX26276628: EndoC-βH1, CVB4-JVB infected, biol rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879324/SRR30879324.lite.1) # SRX26276620: EndoC-βH1, Control, biol rep2, unifected

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879321/SRR30879329.lite.1 # SRX26276629: EndoC-βH1, CVB4-E2 infected, biol rep2

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-7/SRR030/30879/SRR30879326/SRR30879326.lite.1) # SRX26276626: EndoC-βH1, CVB4-JVB infected, biol rep2

fastq-dump --split-files *lite.1 # convert to fastq
