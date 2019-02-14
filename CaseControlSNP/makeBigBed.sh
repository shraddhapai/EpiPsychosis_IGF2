#!/bin/bash

# make a big wig file
BEDSORT=/home/g/gbader/spai/software/kent_utilities/bedSort
BED2BIGBED=/home/g/gbader/spai/software/kent_utilities/bedToBigBed
FETCHCHROMSIZES=/home/g/gbader/spai/software/kent_utilities/fetchChromSizes

snpFile=/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399/hg19/SUS19399.hg19.sorted_-clean.bim
outDir=`pwd`
baseF=`basename $snpFile`

head $snpFile

oF=${outDir}/${baseF}.bed
cat $snpFile | awk '{print "chr"$1"\t"$4"\t"$4"\t"$2}' > $oF
$BEDSORT $oF $oF

#$FETCHCHROMSIZES hg19 > hg19.chrom.sizes

$BED2BIGBED $oF hg19.chrom.sizes ${outDir}/${baseF}.bb

