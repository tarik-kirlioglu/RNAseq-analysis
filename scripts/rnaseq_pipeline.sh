#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "Step1:Prefetch SRA Data"

prefetch -O ${reads}/ ${SRA_ACCESSION}

echo "Step2:Extract fastq Files"

for sra in ${reads}/${SRA_ACCESSION}/*/*.sra
do
	fasterq-dump ${sra} -O ${reads}/
done

echo "Step3:FastQC"

fastqc ${reads}/*.fastq -o ${fastqc}/

echo "Step4:MultiQC"

multiqc ${fastqc}/*_fastqc* -o ${multiqc}/

echo "Step5:Trimming"

for infile in ${reads}/*_1.fastq
do
	base=$(basename ${infile} _1.fastq)
	java -jar ${TRIMMOMATIC_JAR} PE \
				${infile} ${reads}/${base}_2.fastq \
                ${base}_1.trim.fastq ${base}_1.untrim.fastq \
                ${base}_2.trim.fastq ${base}_2.untrim.fastq \
                SLIDINGWINDOW:4:20 MINLEN:25 \
                ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10

done

#moving trimmed files

mv *.trim.* ${trimmed}/

echo "Step6:STAR Mapping"

for infile in ${trimmed}/*_1.trim.fastq
do
   base=$(basename ${infile} _1.trim.fastq)
   STAR --genomeDir ${genome_index}/ \
	--runThreadN ${THREADS} \
	--readFilesIn ${infile} ${trimmed}/${base}_2.trim.fastq \
	--outFileNamePrefix ${mapped}/${base} \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outSAMattributes Standard
done

echo "Step7:Calculating counts with featureCounts"

featureCounts -a ${genome_index}/${GTF_FILE} -o gene_counts.out -T ${THREADS} -p --countReadPairs ${mapped}/*.bam
