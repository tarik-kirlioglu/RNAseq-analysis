
# RNAseq Analysis Project

Welcome to the RNAseq Analysis repository! This project is dedicated to the comprehensive analysis of RNA sequencing data. Below, you'll find all the necessary information to understand, install, and utilize the scripts and data provided.

## Overview
RNA sequencing (RNAseq) is a powerful technique for profiling gene expression. This project includes a complete pipeline for processing RNAseq data, from quality control to differential expression analysis. The dataset used in this project is RNAseq data from [GSE165308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165308), sourced from the NCBI Gene Expression Omnibus (GEO) database. The dataset includes a total of 6 samples, representing both GFI1 knockout and control conditions. All raw reads were processed using quality control tools such as FastQC, and alignment was performed using STAR, followed by gene quantification using featureCounts. Differential gene expression analysis was carried out using the `DESeq2` package in R, identifying genes that were significantly differentially expressed between GFI1 knockout and control samples. Gene Ontology (GO) and KEGG pathway enrichment analyses were also performed using the `clusterProfiler` package in R, allowing us to interpret the biological significance of the differentially expressed genes.

## Getting Started
### Prerequisites
### [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
FastQC is used for quality control of raw sequencing data. You can download and install FastQC from here.

```bash
sudo apt-get install fastqc
```
### [MultiQC](https://github.com/MultiQC/MultiQC)
MultiQC aggregates its results into a single report.You can download Trimmomatic from here.
```bash
sudo apt-get install multiqc
```
### [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
Trimmomatic is used for trimming and cleaning sequencing reads. You can download Trimmomatic from here.
```bash
sudo apt-get install trimmomatic
```
### [SRA Toolkit](https://github.com/ncbi/sra-tools)
SRA Toolkit is used for accessing data from the NCBI Sequence Read Archive (SRA). You can download and install SRA Toolkit from here.

```bash
sudo apt-get install sra-toolkit
```
### [STAR](https://github.com/alexdobin/STAR)
STAR is used for aligning RNAseq reads to a reference genome. You can download and install STAR from here.

```bash
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
make
```
### [featureCounts](https://subread.sourceforge.net/featureCounts.html)
featureCounts is used for quantifying the number of reads mapped to genomic features. It is part of the Subread package, which can be downloaded and installed from here.
```bash
sudo apt-get install subread
```
Ensure you have R installed along with the required packages:
```R
install.packages(c("tidyverse", "ggplot2", "pheatmap"))
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db"))
```
### Installation
Clone the repository to your local machine:
```bash
git clone https://github.com/tarik-kirlioglu/RNAseq-analysis.git
cd RNAseq-analysis
```
## Analysis Workflow
1. **Prefetch Datasets**: Downloading datasets by sra-toolkit.
2. **Quality Control**: Ensuring data quality using tools like FastQC and aggregating the results with MultiQC.
3. **Trimming**: Trimming and cleaning sequencing reads with Trimmomatic
4. **Alignment**: Mapping reads to a reference genome using STAR.
5. **Quantification**: Counting reads with featureCounts.
6. **Statistical Analysis**: Analyzing gene expression with DESeq2.

### Contributions
Contributions are welcome! Feel free to fork this repository and submit pull requests. For major changes, please open an issue first to discuss your proposed changes.

### Contact
If you have any questions or suggestions please feel free to contact me.
