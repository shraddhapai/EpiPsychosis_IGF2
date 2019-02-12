#!/bin/bash

inDir=/home/shraddhapai/Epigenetics/NARSAD/input_files/SeqCapEpi2/LABV_20180618_Capture/trimmed

outF=trimmed_countLines.log
cat /dev/null > $outF
for f in ${inDir}/*log;do
	x=`grep "Both Surviving" $f | awk '{print $7}'`
	baseF=`basename $f`
	echo $baseF

	echo -e "$baseF\t$x" >> $outF
done

