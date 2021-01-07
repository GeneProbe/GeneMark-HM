## RealMetagenomes

Four metagenome assemblies were downloaded from GenBank. All contigs shorter than 1,500 nt were filtered out. 

Skin metagenome
```
mkdir skin
cd skin
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/297/495/GCA_013297495.1_ASM1329749v1/GCA_013297495.1_ASM1329749v1_genomic.fna.gz
gunzip GCA_013297495.1_ASM1329749v1_genomic.fna.gz
probuild --reformat_fasta --in GCA_013297495.1_ASM1329749v1_genomic.fna --out meta.fna  --min_contig 1500 --uppercase 1 --letters_per_line 60 --original
```
Gut metagenome
```
mkdir gut_small
cd gut_small
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/415/335/GCA_900415335.1_SRR3582151/GCA_900415335.1_SRR3582151_genomic.fna.gz
gunzip GCA_900415335.1_SRR3582151_genomic.fna.gz
probuild --reformat_fasta --in GCA_900415335.1_SRR3582151_genomic.fna --out meta.fna --min_contig 1500 --uppercase 1 --letters_per_line 60 --original
```
Gut metagenome
```
mkdir gut_large
cd gut_large
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/604/395/GCA_003604395.1_ASM360439v1/GCA_003604395.1_ASM360439v1_genomic.fna.gz
gunzip GCA_003604395.1_ASM360439v1_genomic.fna.gz
probuild --reformat_fasta --in GCA_003604395.1_ASM360439v1_genomic.fna --out meta.fna --min_contig 1500 --uppercase 1 --letters_per_line 60 --original
```
Oral metagenome
```
mkdir oral
cd oral
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/640/265/GCA_003640265.1_JCVI_supra_1.0/GCA_003640265.1_JCVI_supra_1.0_genomic.fna.gz
gunzip GCA_003640265.1_JCVI_supra_1.0_genomic.fna.gz
probuild --reformat_fasta --in GCA_003640265.1_JCVI_supra_1.0_genomic.fna --out meta.fna --min_contig 1500 --uppercase 1 --letters_per_line 60 --original
```


