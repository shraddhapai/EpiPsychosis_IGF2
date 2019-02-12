#!/bin/bash

alignDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/clean_bam
cd $alignDir
for f in 23Redoneg*.CLEAN.bam;do
	echo $f
samtools flagstat $f > ${f}.stats
done

outF=aligned_CLEAN_reads.log
cat /dev/null > $outF
echo $outF
for f in ${alignDir}/*.CLEAN.bam.stats;do
	x=`grep "properly paired" $f | awk '{print $1}'`
	baseF=`basename $f .CLEAN.bam.stats`
	echo $baseF
	echo $x
	echo -e "$baseF\t$x" >> $outF
done

###for f in *.rmdups.bam;do
###	echo $f
###	samtools flagstat $f > ${f}.stats
###done
###
###for f in *.rmdups.filter.bam;do
###	echo $f
###	samtools flagstat $f > ${f}.stats
###done
