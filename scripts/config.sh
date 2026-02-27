#!/bin/bash

# ============================================================
# Configuration file for RNAseq analysis pipeline
# Edit the paths below according to your system before running.
# ============================================================

# Base working directory
BASE_DIR="/mnt/c/users/honor/desktop/homo_sapiens"

# Directory paths
counts="${BASE_DIR}/counts"
fastqc="${BASE_DIR}/fastqc"
genome_index="${BASE_DIR}/genome_index"
mapped="${BASE_DIR}/mapped"
reads="${BASE_DIR}/reads"
trimmed="${BASE_DIR}/trimmed"
multiqc="${BASE_DIR}/multiqc"

# Trimmomatic
TRIMMOMATIC_JAR="/mnt/c/users/honor/desktop/software/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar"
TRIMMOMATIC_ADAPTERS="/mnt/c/users/honor/desktop/software/Trimmomatic/adapters/TruSeq3-PE-2.fa"

# Reference genome files
GENOME_FASTA="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="Homo_sapiens.GRCh38.112.gtf"

# Thread count
THREADS=8

# SRA project accession
SRA_ACCESSION="PRJNA694054"
