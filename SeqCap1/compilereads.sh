#!/bin/bash

# compile align stats
alignDir=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/output/align

cd $alignDir

###outF=align.RAW.stats.txt
###cat /dev/null > $outF
###for f in *merged.bam.stats;do
###	s=`grep "with itself" $f | cut -f1 -d' '`
###	baseF=`basename $f .merged.bam.stats`
###	echo -e "$baseF\t$s" >> $outF
###done

###outF=align.CLEAN.stats.txt
###cat /dev/null > $outF
###for f in *CLEAN.bam.stats;do
###	s=`grep "with itself" $f | cut -f1 -d' '`
###	baseF=`basename $f .CLEAN.bam.stats`
###	echo -e "$baseF\t$s" >> $outF
###done


outF=align.rmdups.stats.txt
cat /dev/null > $outF
for f in *rmdups.bam.stats;do
	s=`grep "with itself" $f | cut -f1 -d' '`
	baseF=`basename $f .rmdups.bam.stats`
	echo -e "$baseF\t$s" >> $outF
done

outF=align.rmdups.filter.stats.txt
cat /dev/null > $outF
for f in *rmdups.filter.bam.stats;do
	s=`grep "with itself" $f | cut -f1 -d' '`
	baseF=`basename $f .rmdups.filter.bam.stats`
	echo -e "$baseF\t$s" >> $outF
done
