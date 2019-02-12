#!/bin/bash

# align BSseq reads with BSmap
# BSMap params: https://github.com/genome-vendor/bsmap/blob/master/README.txt
BSMAP=/home/shraddhapai/software/bsmap-2.74/bsmap
grch38_fa=/home/shraddhapai/genome_annotation/grch38/fasta/hg38NoAlt_Lambda.fa 
numCores=22; # Dell workstation
rootDir=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi

rawDir=${rootDir}/output/raw
alignDir=${rootDir}/output/align
mkdir -p $alignDir

for r1 in ${rawDir}/test.R1.fastq.gz; do ##*R1*.trim.fastq.gz; do
 	baseF=`basename $r1`
    echo $baseF
    r2Base=${baseF/R1/R2}
    echo $r2Base
	r2=${rawDir}/${r2Base}
	
	#oF=`basename $r1 _R1_001.trim.fastq.gz`
	oF=`basename $r1 _R1_001.trim.fastq.gz`
	oF=${alignDir}/${oF}.bam
	echo "Outfile is $oF"
	# -r 0 : report only unique hits/pairs only
	# -s 16 : seed length
	# -n 1: Map R1 and R2 to both + and - strand (Cokus protocol)
	cmd="$BSMAP -r 0 -s 16 -n 1 -a $r1 -b $r2 -d $grch38_fa -p 24 -o $oF"
	echo $cmd
	`$cmd`
done
