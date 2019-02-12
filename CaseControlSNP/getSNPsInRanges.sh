#!/bin/bash

#' takes bim file and extracts snps in provided genomic ranges
#' to be used to extract SNPs for mQTL analysis.

snpFile=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/plinkQC/postimpute/postimputeqc
dt=171207 #`date +%y%m%d`
outDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/meQTL_postimpute/${dt}

credSNPFile=/home/shraddhapai/Epigenetics/NARSAD/anno/WonGeschwind_2016_Nature_S24.txt
PGC2_snps=/home/shraddhapai/Epigenetics/NARSAD/anno/scz2_plt1e-9.txt
dmpFile=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/out_agesexPMI_PC12_171129.regions-p.bed

BEDTOOLS=/home/shraddhapai/software/bedtools2/bin/bedtools
BEDSORT=/home/shraddhapai/software/kent_utilities/bedSort
PLINK=/home/shraddhapai/software/plink/plink
# ---------------------------------------------------------
mkdir -p $outDir

baseF=`basename $snpFile`

###echo "* Make snp bed file"
snpBed=${outDir}/${baseF}.bed
###cat ${snpFile}.bim | awk '{print "chr"$1"\t"$4"\t"$4"\t"$2}' > $snpBed
###echo "* Sorting"
###echo -e "\tsnps"
###$BEDSORT $snpBed $snpBed
###
#### --------------------------------
#### Now find the snp names in psycharray, that correspond to credible
#### snps and index scz2 snps
###echo "* Resolve names for pgc2 index"
###cat $PGC2_snps | awk '{print $2"\t"$3"\t"$3"\t"$1}' > ${outDir}/scz2_PGC.bed
###$BEDTOOLS intersect -wa -wb -a $snpBed -b ${outDir}/scz2_PGC.bed > ${outDir}/PsychArray_PGC2.bed
###cat ${outDir}/PsychArray_PGC2.bed | cut -f 4 | sort -k1 | uniq > ${outDir}/PsychArray_PGC2.snps
###
###exit 0

###echo "* now do same for credible snps"
###cat $credSNPFile | awk '(NR>1){print $1"\t"$2"\t"$2"\t"$3}' > ${outDir}/scz2_credible.bed
###$BEDTOOLS intersect -wa -wb -a $snpBed -b ${outDir}/scz2_credible.bed > ${outDir}/PsychArray_credible.bed
###cat ${outDir}/PsychArray_credible.bed | cut -f 4 | sort -k1 | uniq > ${outDir}/PsychArray_credible.snps
###
###exit 0 
###echo "* Extend DMP to 1Mb window"
###Rscript extendDMP.R $dmpFile
###dmpDir=`dirname $dmpFile`
###dmpZones=${dmpDir}/dmpQ0.05.1Mb.win.txt
###
###echo "* now do same for snps in dmp zones"
###$BEDSORT $dmpZones $dmpZones
###$BEDTOOLS intersect -wa -wb -a $snpBed -b $dmpZones > ${outDir}/PsychArray_dmp.bed
###cat ${outDir}/PsychArray_dmp.bed | cut -f 4 | sort -k1 | uniq > ${outDir}/PsychArray_dmp.snps
###
###
###echo "* combining (dmp + pgc2 + credible) snps"
###cat ${outDir}/PsychArray_dmp.snps ${outDir}/PsychArray_credible.snps ${outDir}/PsychArray_PGC2.snps | sort -k1 | uniq > ${outDir}/snps2get


#echo "* getting genotypes"
#$PLINK --bfile $snpFile --extract ${outDir}/snps2get \
#    --recode A-transpose --out ${outDir}/geno.meQTL.snps.${dt}

echo "* Listing snps not in LD"
befFile=${outDir}/geno.meQTL.snps.unpruned.${dt}
aftFile=${outDir}/geno.meQTL.snps.pruned.${dt}
$PLINK --bfile $snpFile --extract ${outDir}/snps2get \
  --make-bed --out $befFile
$PLINK --bfile $befFile --exclude high-LD-regions.txt --indep-pairwise 50 5 0.2 --make-bed --out ${outDir}/prune
$PLINK --bfile ${outDir}/geno.meQTL.snps.unpruned.${dt} --extract ${outDir}/prune.prune.in --recode A-transpose --out $aftFile

###echo "* Writing snp source"
###oF=${outDir}/geno.meQTL.snpinfo.${dt}
###cat /dev/null > $oF;
###cat ${outDir}/PsychArray_PGC2.snps | awk '{print $1"\tPGC2 index snp (p<1e-9)"}' >> $oF
###cat ${outDir}/PsychArray_credible.snps | awk '{print $1"\tWonGeschwind SCZ2 credible snps"}' >> $oF
###cat ${outDir}/PsychArray_dmp.snps | awk '{print $1"\tSNPs in DMP window"}' >> $oF
###
