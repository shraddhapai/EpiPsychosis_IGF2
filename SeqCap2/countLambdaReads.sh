#!/bin/bash
# get BS non-conversion rate in all samples
# All C's in lambda phage DNA must be converted to Ts, as these should all be
#' unmethylated.
# Therefore non-conversion rate is (C/CT) and conversion rate is 1-(C/CT).

alignDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/methylation
outFile=${alignDir}/Lambda_readcount.txt

cd $alignDir
cat /dev/null > $outFile

for inFile in *methylation_results.txt.gz;do
	baseF=`basename $inFile .methylation_results.txt.gz`
	echo $baseF
	lc=`zcat $inFile | grep "Lambda" | wc -l`
	echo -e "$baseF\t$lc" >> $outFile
done
