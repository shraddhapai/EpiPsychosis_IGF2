#!/bin/bash

#vcfFile=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlRNAseq/var_calls/BIPOL_02_.vcf
vcfDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlRNAseq/var_calls

PLINK=/home/shraddhapai/software/plink/plink
genoFile=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/final_hg38/SUS19399.hg19.sorted_-clean.hg38

###$PLINK --vcf $vcfFile --make-bed --out ${vcfFile}.plink.full

# run geno2rna.R to get alleles that match between rnaseq vcf and genotype
#$PLINK --bfile $genoFile --make-bed --out tmp
#$PLINK --bfile tmp --extract range set-ranges_geno.txt --recode A-transpose --out ${genoFile}.rnamatch

for  cur in ${vcfDir}/*.vcf; do
	# create set-range
	echo `basename $cur`
	$PLINK --vcf $cur  --extract range set-ranges_rna.txt --recode A-transpose --out ${cur}.plink
done


