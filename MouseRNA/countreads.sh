#!/bin/bash

# count mapped reads in bam file
alDir=/home/shraddhapai/Epigenetics/NARSAD/input_files/MouseRNAseq2/LABV_20180620_RNA/star

outF=${alDir}/readstats.txt
cat /dev/null > $outF
for f in ${alDir}/2[1234567]*_Log.final.out;do
	baseF=`basename $f`
	echo $baseF
	#readcount=`samtools flagstat $f | grep "0 mapped"`
	readcount=`cat $f | grep "Uniquely mapped reads number"`
	echo $readcount
	echo -e "$baseF\t$readcount" >> $outF
done
