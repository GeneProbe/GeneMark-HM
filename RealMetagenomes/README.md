## RealMetagenomes

Four metagenome assemblies were downloaded from GenBank. All contigs shorter than 1,500 nt were filtered out. 

Skin metagenome
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/297/495/GCA_013297495.1_ASM1329749v1/GCA_013297495.1_ASM1329749v1_genomic.fna.gz
gunzip GCA_013297495.1_ASM1329749v1_genomic.fna.gz
mkdir skin
mv GCA_013297495.1_ASM1329749v1_genomic.fna skin
probuild --reformat_fasta --in skin/GCA_013297495.1_ASM1329749v1_genomic.fna --out skin/meta.fna  --min_contig 1500 --uppercase 1 --letters_per_line 60  --original
```
Gut metagenome
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/415/335/GCA_900415335.1_SRR3582151/GCA_900415335.1_SRR3582151_genomic.fna.gz
gunzip GCA_900415335.1_SRR3582151_genomic.fna.gz
mkdir gut_small
mv GCA_900415335.1_SRR3582151_genomic.fna gut_small
probuild --reformat_fasta --in gut_small/GCA_900415335.1_SRR3582151_genomic.fna --out gut_small/meta.fna  --min_contig 1500 --uppercase 1 --letters_per_line 60  --original
```
Gut metagenome
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/604/395/GCA_003604395.1_ASM360439v1/GCA_003604395.1_ASM360439v1_genomic.fna.gz
gunzip GCA_003604395.1_ASM360439v1_genomic.fna.gz
mkdir gut_large
mv GCA_003604395.1_ASM360439v1_genomic.fna gut_large
probuild --reformat_fasta --in gut_large/GCA_003604395.1_ASM360439v1_genomic.fna --out gut_large/meta.fna  --min_contig 1500 --uppercase 1 --letters_per_line 60  --original
```
Oral metagenome
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/640/265/GCA_003640265.1_JCVI_supra_1.0/GCA_003640265.1_JCVI_supra_1.0_genomic.fna.gz
gunzip GCA_003640265.1_JCVI_supra_1.0_genomic.fna.gz
mkdir oral
mv GCA_003640265.1_JCVI_supra_1.0_genomic.fna oral
probuild --reformat_fasta --in oral/GCA_003640265.1_JCVI_supra_1.0_genomic.fna --out oral/meta.fna  --min_contig 1500 --uppercase 1 --letters_per_line 60  --original
```


