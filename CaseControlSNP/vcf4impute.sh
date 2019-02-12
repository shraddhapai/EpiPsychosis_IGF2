#!/bin/bash

# convert snp data to vcf format for imputation on Michigan imputation server
inFile=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/plinkQC/hg19/SUS19399.hg19
PLINK=/home/shraddhapai/software/plink/plink
BGZIP=/home/shraddhapai/software/htslib-1.5/bgzip

chrSet=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22);

outDir=`dirname $inFile`/impute
mkdir -p $outDir
baseF=`basename $inFile`
outF=${outDir}/${baseF}
for chrom in "${chrSet[@]}"; do
	echo $chrom
	$PLINK --bfile $inFile --chr $chrom --recode vcf --out ${outF}_chr${chrom}
	vcf-sort ${outF}_chr${chrom}.vcf | $BGZIP -c > ${outF}_chr${chrom}.vcf.gz
	
	
done
