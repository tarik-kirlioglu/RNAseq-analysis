#!/bin/bash

counts="/mnt/c/users/honor/desktop/c_elegans/counts/"
fastqc="/mnt/c/users/honor/desktop/c_elegans/fastqc/"
genome_index="/mnt/c/users/honor/desktop/c_elegans/genome_index/"
mapped="/mnt/c/users/honor/desktop/c_elegans/mapped/"
reads="/mnt/c/users/honor/desktop/c_elegans/reads/"
trimmed="/mnt/c/users/honor/desktop/c_elegans/trimmed/"

echo "Step1:Prefetch SRA Data"
if false 
then
prefetch -O ${reads} SRR6822884 SRR6822885

echo "Step2:Extract fastq Files"

fasterq-dump ${reads}*.sra -O ${reads}


echo "Step3:FastQC"

fastqc ${reads}*.fastq -o ${fastqc}


echo "Step4:Trimming"

for infile in ${reads}*_1.fastq
do
	base=$(basename ${infile} _1.fastq)
	java -jar /mnt/c/users/honor/desktop/software/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE \
				${infile} ${reads}${base}_2.fastq \
                ${base}_1.trim.fastq ${base}_1.untrim.fastq \
                ${base}_2.trim.fastq ${base}_2.untrim.fastq \
                SLIDINGWINDOW:4:20 MINLEN:25 \
                ILLUMINACLIP:/mnt/c/users/honor/desktop/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 

done


#moving trimmed files

mv *.trim.* ${trimmed}

#downloading gtf and genom fasta files

wget https://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.111.gtf.gz

echo "Step5:STAR Genome Index"

STAR --runMode genomeGenerate \
	 --genomeDir ${genome_index} \
	 --genomeFastaFiles Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
	 --sjdbGTFfile Caenorhabditis_elegans.WBcel235.111.gtf \
	 --runThreadN 2



echo "Step6:STAR Mapping"


for infile in ${trimmed}*_1.trim.fastq
do
   base=$(basename ${infile} _1.trim.fastq)
   STAR --genomeDir ${genome_index} \
	--runThreadN 2 \
	--readFilesIn ${infile} ${trimmed}${base}_2.trim.fastq \
	--outFileNamePrefix ${mapped}${base} \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outSAMattributes Standard
done

fi

featureCounts -a Caenorhabditis_elegans.WBcel235.111.gtf  -o count.out -T 2 -p  ${mapped}*.bam

