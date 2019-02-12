#!/bin/bash

METHRATIO=/home/shraddhapai/software/bsmap-2.74/methratio.py
#PICARD=/home/shraddhapai/software/picard.jar
PICARD=/home/shraddhapai/software/picard-tools/1.98/picard-1.98.jar
SAMTOOLS=/home/shraddhapai/software/samtools-1.5/samtools
BAMUTIL=/home/shraddhapai/software/bamUtil/bin/bam
BEDTOOLS=/home/shraddhapai/software/bedtools2/bin/bedtools
BGZIP=/home/shraddhapai/software/htslib-1.5/bgzip
TABIX=/home/shraddhapai/software/htslib-1.5/tabix
hg38_fa=/home/shraddhapai/genome_annotation/grch38/fasta/hg38NoAlt_Lambda.fa
hg38_chromSizes=/home/shraddhapai/genome_annotation/grch38/hg38.chrom.sizes
primaryTargets=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/input/OID43737_hg19_161128_design_deliverables/OID43737_hg19_161128_primary_targets.hg38.bed
captureTargets=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/input/OID43737_hg19_161128_design_deliverables/OID43737_hg19_161128_capture_targets.hg38.bed

sampName=$1;
alignDir=$2;

cd $alignDir
baseF=`basename $sampName .bam`
### ----------------------------------------------------
### Combine reads sequenced in different lanes
### ----------------------------------------------------
echo "Started post-alignment processing"
date

###	echo "Merging"
###	fNames=`ls $sampName*L00*.bam`
###	echo $fNames
###	cmd="$SAMTOOLS merge ${sampName}.merged.bam $fNames"
###	echo `$cmd`

### ----------------------------------------------------
### Remove duplicates for top/bottom strands separately
### and then remerge.
### ----------------------------------------------------

	echo "split by top/bottom"
 
	$SAMTOOLS view -H ${sampName}.merged.bam > ${baseF}.header.txt
	$SAMTOOLS view -S ${sampName}.merged.bam | grep "ZS:Z:+" > ${baseF}.top.sam.tmp
 cat ${baseF}.header.txt ${baseF}.top.sam.tmp > ${baseF}.top.sam

	$SAMTOOLS view -S ${sampName}.merged.bam | grep "ZS:Z:-" > ${baseF}.bottom.sam.tmp
 cat ${baseF}.header.txt ${baseF}.bottom.sam.tmp > ${baseF}.bottom.sam
	
	rm ${baseF}.top.sam.tmp
	rm ${baseF}.bottom.sam.tmp
	rm ${baseF}.header.txt

	# now sort separately
	echo "sort each"
	$SAMTOOLS sort ${baseF}.top.sam > ${baseF}.top.sorted.bam
	$SAMTOOLS sort ${baseF}.bottom.sam > ${baseF}.bottom.sorted.bam

	# clean sam files
	rm ${baseF}.top.sam ${baseF}.bottom.sam
	
	# remove duplicates
	echo "remove duplicates for each"
	java -Xmx4g -Xms4g -jar $PICARD MarkDuplicates \
		VALIDATION_STRINGENCY=LENIENT \
		INPUT=${baseF}.top.sorted.bam \
		OUTPUT=${baseF}.top.rmdups.bam \
		METRICS_FILE=${baseF}.top.rmdups_metrics.txt \
		REMOVE_DUPLICATES=true ASSUME_SORTED=true \
		CREATE_INDEX=true

	java -Xmx4g -Xms4g -jar $PICARD MarkDuplicates \
		VALIDATION_STRINGENCY=LENIENT \
		INPUT=${baseF}.bottom.sorted.bam \
		OUTPUT=${baseF}.bottom.rmdups.bam \
		METRICS_FILE=${baseF}.bottom.rmdups_metrics.txt \
		REMOVE_DUPLICATES=true ASSUME_SORTED=true \
		CREATE_INDEX=true

echo "now remerge top/bottom"
$SAMTOOLS merge ${baseF}.rmdups.bam ${baseF}.top.rmdups.bam ${baseF}.bottom.rmdups.bam 

#cleanup again
rm ${baseF}.top.rmdups.bam ${baseF}.bottom.rmdups.bam 
rm ${baseF}.top.sorted.bam ${baseF}.bottom.sorted.bam

### ----------------------------------------------------
### Clean BAM file
### ----------------------------------------------------
	echo "Filter bam file"
	echo "filtering"
	$SAMTOOLS view -F 4 -f 2 -b ${baseF}.rmdups.bam > ${baseF}.rmdups.filter.bam

	echo "clipping overhanging reads"
	$BAMUTIL clipOverlap --stats --in ${baseF}.rmdups.filter.bam --out ${baseF}.CLEAN.bam

	echo "sort and indexing BAM"
	$SAMTOOLS sort ${baseF}.CLEAN.bam ${baseF}.tmp
	mv ${baseF}.tmp.bam ${baseF}.CLEAN.bam
	$SAMTOOLS index ${baseF}.CLEAN.bam

	echo "mapping metrics from picard"
	java -Xmx4g -Xms4g -jar $PICARD \
		CollectAlignmentSummaryMetrics \
		METRIC_ACCUMULATION_LEVEL=ALL_READS \
		INPUT=${baseF}.CLEAN.bam \
		OUTPUT=${baseF}.CLEAN.picard_alignment_metrics.txt \
		REFERENCE_SEQUENCE=$hg38_fa \
		VALIDATION_STRINGENCY=LENIENT

### ----------------------------------------------------
### Compute targeted capture stats
### ----------------------------------------------------

echo "Creating intervals and computing stats"
echo "primary"
cat $primaryTargets | awk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > DESIGN_target_body.txt
echo "capture"
cat $captureTargets  | awk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > DESIGN_bait_body.txt

echo "catting header"
$SAMTOOLS view -H ${baseF}.CLEAN.bam > ${baseF}.header.txt
cat ${baseF}.header.txt DESIGN_target_body.txt > DESIGN_target_intervals.txt
cat ${baseF}.header.txt DESIGN_bait_body.txt > DESIGN_bait_intervals.txt

rm ${baseF}.header.txt

echo "picard calchsmetrics"
java -Xmx4g -Xms4g -jar $PICARD CollectHsMetrics \
	INPUT=${baseF}.CLEAN.bam \
	OUTPUT=${baseF}.picard_hs_metrics.txt \
	R=$hg38_fa \
	BAIT_INTERVALS=DESIGN_bait_intervals.txt \
	TARGET_INTERVALS=DESIGN_target_intervals.txt \
	METRIC_ACCUMULATION_LEVEL=ALL_READS  \
	VALIDATION_STRINGENCY=LENIENT

 echo "8. estimate insert size distribution"
	java -Xmx4g -jar $PICARD CollectInsertSizeMetrics \
		INPUT=${baseF}.CLEAN.bam \
		OUTPUT=${baseF}_picard_insert_size_metrics.txt \
	 	VALIDATION_STRINGENCY=LENIENT \
		HISTOGRAM_FILE=${baseF}_picard_insert_size_plot.pdf

echo "9. add target padding"
$BEDTOOLS slop -i $captureTargets -g $hg38_chromSizes -b 100 | $BEDTOOLS sort -i - | $BEDTOOLS merge -i - > DESIGN_padded_capture_target.bed

####### this step not need
####### 10. sum total region size in target bed
########$BEDTOOLS genomecov -i DESIGN.bed -g chromosome_sizes.txt -max 1 | grep -P "genome\t1" | cut -f 3

###echo "11. get primary target reads"
###$BEDTOOLS intersect -bed -abam ${baseF}.CLEAN.bam -b $primaryTargets > ${baseF}.primary.target.reads.bed
###echo "12. get capture target reads (excluding padding)"
###$BEDTOOLS intersect -bed -abam ${baseF}.CLEAN.bam -b $captureTargets > ${baseF}.capture.target.reads.bed

echo "12. calc depth coverage"
$BEDTOOLS coverage -a $primaryTargets -b ${baseF}.CLEAN.bam > ${baseF}.primary.target.coverage.bed 

### ----------------------------------------------------
### Get methylation
### ----------------------------------------------------
echo "13. Calculate %m"
samDir=`dirname $SAMTOOLS`
python $METHRATIO -d $hg38_fa -s $samDir -m 1 -z -i skip -o ${baseF}.methylation_results.txt ${baseF}.CLEAN.bam

echo "convert to tabix"
$BGZIP ${baseF}.methylation_results.txt 
$TABIX -s 1 -b 2 -e 2 -S 1 ${baseF}.methylation_results.txt.gz

echo "Ended post-alignment processing"
date

### Add later
###echo "calculate bs efficiency"
###samDir=`dirname $SAMTOOLS`
###python $METHRATIO -d $hg38_fa -s $samDir -m 1 -z -i skip -o -c Lambda_J02459.1 ${baseF}.Lambda.methylation_results.txt ${baseF}.CLEAN.bam

