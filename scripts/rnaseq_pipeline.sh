#!/bin/bash

counts="/mnt/c/users/honor/desktop/homo_sapiens/counts"
fastqc="/mnt/c/users/honor/desktop/homo_sapiens/fastqc"
genome_index="/mnt/c/users/honor/desktop/homo_sapiens/genome_index"
mapped="/mnt/c/users/honor/desktop/homo_sapiens/mapped"
reads="/mnt/c/users/honor/desktop/homo_sapiens/reads"
trimmed="/mnt/c/users/honor/desktop/homo_sapiens/trimmed"
multiqc="/mnt/c/users/honor/desktop/homo_sapiens/multiqc"

echo "Step1:Prefetch SRA Data"

prefetch -O ${reads}/ PRJNA694054

echo "Step2:Extract fastq Files"

fasterq-dump ${reads}/*.sra -O ${reads}/

echo "Step3:FastQC"

fastqc ${reads}/*.fastq -o ${fastqc}/

echo "Step4:MultiQC"

multiqc ${fastqc}/*_fastqc* -o ${multiqc}/

echo "Step5:Trimming"

for infile in ${reads}/*_1.fastq
do
	base=$(basename ${infile} _1.fastq)
	java -jar /mnt/c/users/honor/desktop/software/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE \
				${infile} ${reads}/${base}_2.fastq \
                ${base}_1.trim.fastq ${base}_1.untrim.fastq \
                ${base}_2.trim.fastq ${base}_2.untrim.fastq \
                SLIDINGWINDOW:4:20 MINLEN:25 \
                ILLUMINACLIP:/mnt/c/users/honor/desktop/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 

done

#moving trimmed files

mv *.trim.* ${trimmed}/

echo "Step6:STAR Mapping"

for infile in ${trimmed}/*_1.trim.fastq
do
   base=$(basename ${infile} _1.trim.fastq)
   STAR --genomeDir ${genome_index}/ \
	--runThreadN 8 \
	--readFilesIn ${infile} ${trimmed}/${base}_2.trim.fastq \
	--outFileNamePrefix ${mapped}/${base} \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outSAMattributes Standard
done

echo "Step7:Calculating counts with featureCounts"

featureCounts -a Homo_sapiens.GRCh38.112.gtf  -o gene_counts.out -T 8 -p  ${mapped}/*.bam
