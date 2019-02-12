#!/bin/bash

indir=/home/shraddhapai/Epigenetics/NARSAD/input_files/SeqCapEpi2/LABV_20180618_Capture/fastq

outfile=fastq_countLines.log
cat /dev/null >>$outfile
for f in ${indir}/*.fastq.gz;do
	baseF=`basename $f`
	echo $baseF
	nl=`zcat $f | wc -l | awk '{rd=$1 / 4; print rd}'`
	echo -e "$baseF\t$nl" >> $outfile
done
