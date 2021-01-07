# GeneMark-HM
_GeneMark-HM: Improving Gene Prediction in DNA Sequences of Human Microbiome_  
(manuscript submitted for publication, 2020)  

## Overview
GeneMark-HM is a pipeline for protein coding gene prediction in assembled metagenomes from human microbiome (HM). Pipeline integrates three algorithms from GeneMark gene finder family, MetaGeneMark-2, GeneMarkS-2 and Genemark.hmm-2, into a workflow with parameters optimized for analysis of human metagenomes. 

![diagrmm](./docs/diagramm.jpg)

This pipeline is built around a fact, that gene prediction with species specific parameters usually outperforms gene prediction with more general metagenomics parameters. Species specific parameters are impossible to estimate from short contigs. For many species present in human gut microbiome, complete or almost genomes are available in public repositories. Accurate estimate of parameters can be done for such species by GeneMarkS-2 algorithm using complete genomes. This pipeline leverages the availability of data for many taxa. Database of species specific parameters was built for more than five thousand taxa, which are known to be present in the human gut metagenomes. Pipeline checks if a contig can be assigned to one of the reference taxa and in case of reliable assignment uses corresponding taxa specific parameters for gene finding in a contig. In cases when no taxa assignment was made one of the two algorithms GeneMarkS-2 or MetaGeneMark-2 are used to generate the gene prediction. Selection of algorithms on the latest step is based on the contig length.

Two databases form the core part of this pipeline. Database of genes is used on taxa assignment step and database of species specific parameters on gene prediction step. Quality of these two databases determine the gain in gene annotation accuracy from initial MataGeneMark-2 level.

## Pipeline code
* GeneMark-HM software is located in "bin" folder
* pan-genome sequences are located on Amazon S3 storage
* pan-genome specific gene finding models are located on Amazon S3 storage

## Instalation
Installation instructions are available in [INSTALL](INSTALL) file.

## Data and code used in the paper

Table 6. Collection of statistics on genes predicted in real metagenomes [RealMetagenomes](RealMetagenomes). 


