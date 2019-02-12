#!/bin/bash

alignDir=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/output/align
cd $alignDir
###for f in *.merged.bam;do
###	echo $f
###	samtools flagstat $f > ${f}.stats
###done
###for f in *.CLEAN.bam;do
###	echo $f
###	samtools flagstat $f > ${f}.stats
###done

for f in *.rmdups.bam;do
	echo $f
	samtools flagstat $f > ${f}.stats
done

for f in *.rmdups.filter.bam;do
	echo $f
	samtools flagstat $f > ${f}.stats
done
