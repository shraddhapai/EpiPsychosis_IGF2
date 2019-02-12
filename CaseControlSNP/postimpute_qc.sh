#!/bin/bash

# qc from http://www.well.ox.ac.uk/~wrayner/tools/Post-Imputation.html

vcfparse=/home/shraddhapai/software/vcfparse.pl
ic=/home/shraddhapai/software/IC/ic.pl
KGP_REF=/home/shraddhapai/software/preimpute_check/1000GP_Phase3_combined.legend
PLINK=/home/shraddhapai/software/plink/plink


impDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/impute

###password=`cat ${impDir}/password.txt`

###cat /dev/null > ${impDir}/allPlinkFiles.txt
###for chr in {2..22}; do
###	#perl $vcfparse -d $impDir -o ${impDir}/vcfparse -g
###	#perl $ic -d ${impDir}/vcfparse -r $KGP_REF -g -o ${impDir}/vcfparse
###	echo "chrom $chr"
###	
###	echo "unzipping"
###	unzip -P $password -o ${impDir}/chr_$chr.zip -d $impDir
###
###	# filter by info score
###	zcat ${impDir}/chr${chr}.info.gz | awk '(NR>1) {if ($7>=0.7){print}}'| cut -f 1 >  ${impDir}/chr${chr}.info0.7.txt
###
###	echo "Extract high-conf snps" 
###	$PLINK --vcf ${impDir}/chr${chr}.dose.vcf.gz --make-bed --out ${impDir}/chr${chr}
###	$PLINK --bfile ${impDir}/chr${chr} --extract ${impDir}/chr${chr}.info0.7.txt --make-bed --out ${impDir}/chr${chr}.info0.7
###
###	# cleanup
###	rm ${impDir}/chr${chr}.bim ${impDir}/chr${chr}.bed ${impDir}/chr${chr}.fam
###	rm ${impDir}/chr${chr}.info.gz
###	rm ${impDir}/chr${chr}.dose.vcf.gz
###	rm ${impDir}/chr${chr}.dose.vcf.gz.tbi
###
###	echo "${impDir}/chr${chr}.info0.7.bed ${impDir}/chr${chr}.info0.7.bim ${impDir}/chr${chr}.info0.7.fam" >> ${impDir}/allPlinkFiles.txt
###
###done

# manually edited list -- don't just run this line below
#ls ${impDir}/chr*bed > ${impDir}/allPlinkFiles.txt

# -------
# merge across chroms
#### first attempt will throw error because of alleles with 3+ variants
###$PLINK --bfile ${impDir}/chr1.info0.7 --merge-list ${impDir}/allPlinkFiles.txt --make-bed --out ${impDir}/comb_r1
###
####### remove missnps and duplicate vars
###for chr in {1..22};do
###	$PLINK --bfile ${impDir}/chr${chr}.info0.7 --exclude ${impDir}/comb_r1-merge.missnp --make-bed --out ${impDir}/chr${chr}.nomissnp
###	$PLINK --bfile ${impDir}/chr${chr}.nomissnp --list-duplicate-vars
###	$PLINK --bfile ${impDir}/chr${chr}.nomissnp --exclude plink.dupvar --make-bed --out ${impDir}/chr${chr}.nodups
###done

cat /dev/null > ${impDir}/allPlinkFiles.txt
for chr in {2..22};do
	echo "${impDir}/chr${chr}.nodups.bed ${impDir}/chr${chr}.nodups.bim ${impDir}/chr${chr}.nodups.fam" >> ${impDir}/allPlinkFiles.txt
done

# now try merging again
$PLINK --bfile ${impDir}/chr1.nodups --merge-list ${impDir}/allPlinkFiles.txt --make-bed --out ${impDir}/postimputeqc

# then run plinkqc script with these instead of original

# genetic-epigenetic just extract snps in the region


# loop over each file if needed.


# and now you proceed with plink cleaning merging with HM3 etc etc
#$PLINK --bfile ${impDir}/chr11.info0.7 --geno 0.01 --hwe 0.000001 --make-bed -out ${impDir}/chr11.clean

# finally convert to vcf for ASE
#$PLINK --bfile ${impDir}/chr11.clean --extract range ${impDir}/IGF2_ranges.txt --recode vcf --out ${impDir}/chr11.clean.tovcf






