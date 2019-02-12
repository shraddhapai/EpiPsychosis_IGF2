#!/bin/bash

#' allele-specific expression

GATK=/home/shraddhapai/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
PICARD=/home/shraddhapai/software/picard_2.9.4/picard.jar
SAMTOOLS=/home/shraddhapai/software/samtools-1.5/samtools
hg38_fa=/home/shraddhapai/genome_annotation/grch38/fasta/hg38.Lee.fa
hg38_dict=/home/shraddhapai/genome_annotation/grch38/fasta/hg38.Lee.dict
PICARD=/home/shraddhapai/software/picard_2.9.4/picard.jar
SAMTOOLS=/home/shraddhapai/software/samtools-1.5/samtools
PLINK=/home/shraddhapai/software/plink/plink

genoIn=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/plinkQC/postimpute/SUS19399.posimp.hg19.sorted_-clean
outDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlRNAseq/var_calls
bamDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlRNAseq/LeeData/bam_files
snpFile=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/test.vcf

for f in ${bamDir}/CTRL_86*.bam; do
	baseF=`basename $f .bam`
	echo $baseF
	sampName=`basename $baseF Aligned.sortedByCoord.out`
	newF=${outDir}/${baseF}

	echo -e "\tpicard"
	#java -jar $PICARD AddOrReplaceReadGroups I=$f O=$newF.GP.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=1
	echo -e "\tindex"
	#$SAMTOOLS index ${newF}.GP.bam

	echo -e "\treassign mapq"
	#java -jar -Xmx32g $GATK -T SplitNCigarReads -R $hg38_fa -I ${newF}.GP.bam -o ${newF}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

	echo -e "base quality calibration".
#	java -jar -Xmx32g $GATK -T PrintReads -R $hg38_fa -I ${newF}.split.bam -o ${newF}.bqsr.bam
	echo -e "\tindex"
	#$SAMTOOLS index ${newF}.mapqFixed.bam
	echo -e "\thaplocaller"
	java -jar -Xmx50g $GATK -T HaplotypeCaller -nct 16 -R $hg38_fa -I ${newF}.bqsr.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o ${outDir}/${sampName}.vcf
###	$PLINK --vcf ${outDir}/${sampName}.vcf --make-bed --out ${outDir}/${sampName}.plink
done
