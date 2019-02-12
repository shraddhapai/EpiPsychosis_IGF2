#!/bin/bash

# FastQC processing
#FASTQC=/home/spai/Other_Software/FastQC/fastqc
FASTQC=/home/shraddhapai/software/FastQC/fastqc

#rootDir=/scratch/g/gbader/spai/NARSAD2014/input_files/CaseControlRNAseq
#rootDir=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi
#inDir=${rootDir}/input/Samples
rootDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2
inDir=${rootDir}/trimmed
outDir=${rootDir}/fastqc

#-----------------------------
#module load java/8.0

mkdir -p $outDir

#for f in ${inDir}/*re/Files/*fastq.gz;do
for f in ${inDir}/*trim.fastq.gz;do
	$FASTQC $f --outdir=$outDir --extract
done


