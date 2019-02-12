#!/bin/bash

# converts coordinates from input genome build to hg19 coordinates.
# 1000Genomes ref panel is in hg19 so our chip coordinates need to be in same build
# note that we do not convert chrX,Y or M because we do not analyze these chromosomes

liftOver=~/software/kent_utilities/liftOver
LO_CHAIN=~/software/kent_utilities/hg19ToHg38.over.chain.gz
PLINK=~/software/plink/plink

rootDir=~/Epigenetics/NARSAD/output_files/CaseControlSNP
HG38DIR=${rootDir}/final_hg38
HG19DIR=${rootDir}/final_hg19
LIFTDIR=${rootDir}/liftOver

mkdir -p $LIFTDIR
mkdir -p $HG19DIR
mkdir -p $LIFTDIR

baseF=postimputeqc

###cat ${HG19DIR}/${baseF}.bim | awk '{print "chr"$1"\t"$4"\t"$4+1"\t"$2}' > ${LIFTDIR}/lOver.inBuild

### $liftOver ${LIFTDIR}/lOver.inBuild $LO_CHAIN ${LIFTDIR}/lOver.hg38 ${LIFTDIR}/lOver.unmapped

### ONLY do this if some snps have been unmapped in converting
### compile list of snps that didn't map successfully and 
### therefore need to be excluded
###grep -v ^# ${LIFTDIR}/lOver.unmapped | awk '{print $4}' > ${LIFTDIR}/toremove.txt
### convert to map and ped file for reordering
###$PLINK --bfile ${HG19DIR}/${baseF} --exclude ${LIFTDIR}/toremove.txt --recode --out ${LIFTDIR}/data.inBuild

###exit 0
###
######## compare SNPs in map file with the order in the liftover output. 
######## If the order is the same, we do not need to specially reorder the columns
###cut -f4 ${LIFTDIR}/lOver.hg38 > ${LIFTDIR}/a.snp
###cut -f2 ${LIFTDIR}/data.inBuild.map > ${LIFTDIR}/b.snp
###curwd=`pwd`
###cd $LIFTDIR
###DIFF=$(diff a.snp b.snp)
###cd $curwd
###if [ "$DIFF" != "" ] ; then
###	echo "*** Order of SNPs in LiftOver doesn't match that in plink ***"
###	echo "*** Manual reordering may be required                     ***"
###	exit 0
###else
###	echo "   * SNP order in LiftOver and plink match!"
###	echo "   * Proceed with conversion to bed"
###	rm ${LIFTDIR}/a.snp ${LIFTDIR}/b.snp
###fi
###
####### now that we have confirmed order is correct, convert back to 
###### binary format.
$PLINK --file ${LIFTDIR}/data.inBuild --make-bed --out ${LIFTDIR}/${baseF}.hg38
paste ${LIFTDIR}/${baseF}.hg38.bim ${LIFTDIR}/lOver.hg38 > ${LIFTDIR}/tmp.bim
cat ${LIFTDIR}/tmp.bim | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > ${LIFTDIR}/tmp2.bim
mv ${LIFTDIR}/tmp2.bim ${LIFTDIR}/${baseF}.hg38.bim

rm ${LIFTDIR}/tmp.bim ${LIFTDIR}/data.inBuild.*

