
# RNAseq Analysis Project

Welcome to the RNAseq Analysis repository! This project is dedicated to the comprehensive analysis of RNA sequencing data. Below, you'll find all the necessary information to understand, install, and utilize the scripts and data provided.

## Overview
RNA sequencing (RNAseq) is a powerful technique for profiling gene expression. This project includes a complete pipeline for processing RNAseq data, from quality control to differential expression analysis.

## Getting Started
### Prerequisites
Ensure you have R installed along with the required packages:
```R
install.packages(c("tidyverse", "ggplot2"))
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "pheatmap", "EnhancedVolcano, "org.Hs.eg.db"))
```
### Installation
Clone the repository to your local machine:
```bash
git clone https://github.com/tarik-kirlioglu/RNAseq-analysis.git
cd RNAseq-analysis
```
## Analysis Workflow
1. **Quality Control**: Using tools like FastQC to ensure data quality.
2. **Alignment**: Mapping reads to a reference genome using STAR.
3. **Quantification**: Counting reads with featureCounts.
4. **Statistical Analysis**: Analyzing gene expression with DESeq2.

### Contributions
Contributions are welcome! Feel free to fork this repository and submit pull requests. For major changes, please open an issue first to discuss your proposed changes.

### Contact
If you have any questions or suggestions, please feel free to contact Tarık Kırlıoğlu.
