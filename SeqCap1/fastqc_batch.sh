#!/bin/bash

# FastQC processing
#FASTQC=/home/spai/Other_Software/FastQC/fastqc
FASTQC=/home/shraddhapai/software/FastQC/fastqc

#rootDir=/scratch/g/gbader/spai/NARSAD2014/input_files/CaseControlRNAseq
rootDir=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi
inDir=${rootDir}/input/Samples
outDir=${rootDir}/fastqc

#-----------------------------
#module load java/8.0

mkdir -p $outDir

for f in ${inDir}/*re/Files/*fastq.gz;do
	$FASTQC $f --outdir=$outDir --extract
done


