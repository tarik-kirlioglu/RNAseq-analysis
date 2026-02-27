#!/bin/bash

# Genome Setup Script
# This script downloads the reference genome and GTF annotation files,
# then builds the STAR genome index. Run this once before the main pipeline.

genome_index="/mnt/c/users/honor/desktop/homo_sapiens/genome_index"

echo "Step1: Downloading reference genome and GTF annotation files"

wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz

echo "Step2: Extracting files"

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz

echo "Step3: Building STAR Genome Index"

STAR --runMode genomeGenerate \
	 --genomeDir ${genome_index}/ \
	 --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	 --sjdbGTFfile Homo_sapiens.GRCh38.112.gtf \
	 --runThreadN 8

echo "Genome setup completed!"
