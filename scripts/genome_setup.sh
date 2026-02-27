#!/bin/bash

# Genome Setup Script
# This script downloads the reference genome and GTF annotation files,
# then builds the STAR genome index. Run this once before the main pipeline.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "Step1: Downloading reference genome and GTF annotation files"

wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/${GENOME_FASTA}.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/${GTF_FILE}.gz

echo "Step2: Extracting files"

gunzip ${GENOME_FASTA}.gz
gunzip ${GTF_FILE}.gz

echo "Step3: Building STAR Genome Index"

STAR --runMode genomeGenerate \
	 --genomeDir ${genome_index}/ \
	 --genomeFastaFiles ${GENOME_FASTA} \
	 --sjdbGTFfile ${GTF_FILE} \
	 --runThreadN ${THREADS}

echo "Genome setup completed!"
