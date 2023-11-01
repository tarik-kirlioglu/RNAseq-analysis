#!/bin/bash

ref="/mnt/c/users/honor/desktop/rnaseq_pipeline/ref"
reads="/mnt/c/users/honor/desktop/rnaseq_pipeline/reads"
results="/mnt/c/users/honor/desktop/rnaseq_pipeline/results"
quality="/mnt/c/users/honor/desktop/rnaseq_pipeline/fastqc"
mapStar="/mnt/c/users/honor/desktop/rnaseq_pipeline/mapStar"
trim="/mnt/c/users/honor/desktop/rnaseq_pipeline/trim"

echo "STEP 1:Prefetch SRA Data"

prefetch -O ${reads}/ SRR12439141 SRR12439142 SRR12439143 SRR12439144 SRR12439145 SRR12439146

echo "STEP 2:Fasterq-dump SRA Data"

fasterq-dump ${reads}/*.sra -O ${reads}/

echo "STEP 3:FastQC"

fastqc ${reads}/*.fastq -o ${quality}/

# NOTE: If fastqc gives bad results then it must be trimming 
echo "STEP 4:Trimmomatic"

java -jar /mnt/c/users/honor/desktop/software/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE \
	${reads}/SRR12439141_1.fastq \
	${reads}/SRR12439141_2.fastq \
	${trim}/SRR12439141_1P.fastq ${trim}/SRR12439141_1U.fastq \
	${trim}/SRR12439141_2P.fastq ${trim}/SRR12439141_2U.fastq \
	ILLUMINACLIP:/mnt/c/users/honor/desktop/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
	LEADING:3 TRAILING:3 MINLEN:36

if false 
then
	
echo "STEP 5:STAR Genome Index"

STAR --runMode genomeGenerate \
	--genomeDir ${ref}/ \
	--genomeFastaFiles ${ref}/Pvulgaris_442_v2.1.softmasked.fa \
	--sjdbGTFfile ${ref}/Pvulgaris_442_v2.1.gene_exons.gtf \
	--runThreadN 2 \
	--genomeSAindexNbases 12

echo "STEP 5:STAR Alignment"

STAR --genomeDir ${ref}/ \
	--runThreadN 2 \
	--readFilesIn ${reads}/SRR12439141_1.fastq ${reads}/SRR12439141_2.fastq \
	--outFileNamePrefix ${mapStar}/SRR12439141 \ 
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--quantMode GeneCounts \
	--outSAMattributes Standard

STAR --genomeDir ${ref}/ \
	--runThreadN 2 \
	--readFilesIn ${reads}/SRR12439142_1.fastq ${reads}/SRR12439142_2.fastq \
	--outFileNamePrefix ${mapStar}/SRR12439142 \ 
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--quantMode GeneCounts \
	--outSAMattributes Standard

STAR --genomeDir ${ref}/ \
	--runThreadN 2 \
	--readFilesIn ${reads}/SRR12439143_1.fastq ${reads}/SRR12439143_2.fastq \
	--outFileNamePrefix ${mapStar}/SRR12439143 \ 
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--quantMode GeneCounts \
	--outSAMattributes Standard

STAR --genomeDir ${ref}/ \
	--runThreadN 2 \
	--readFilesIn ${reads}/SRR12439144_1.fastq ${reads}/SRR12439144_2.fastq \
	--outFileNamePrefix ${mapStar}/SRR12439144 \ 
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--quantMode GeneCounts \
	--outSAMattributes Standard

STAR --genomeDir ${ref}/ \
	--runThreadN 2 \
	--readFilesIn ${reads}/SRR12439145_1.fastq ${reads}/SRR12439145_2.fastq \
	--outFileNamePrefix ${mapStar}/SRR12439145 \ 
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--quantMode GeneCounts \
	--outSAMattributes Standard

STAR --genomeDir ${ref}/ \
	--runThreadN 2 \
	--readFilesIn ${reads}/SRR12439146_1.fastq ${reads}/SRR12439146_2.fastq \
	--outFileNamePrefix ${mapStar}/SRR12439146 \ 
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--quantMode GeneCounts \
	--outSAMattributes Standard

echo "STEP 6: featureCounts"

featureCounts -a ${ref}/Pvulgaris_442_v2.1.gene_exons.gtf -o ${results}/count.out -T 2 -p  ${mapStar}/*.bam

fi

