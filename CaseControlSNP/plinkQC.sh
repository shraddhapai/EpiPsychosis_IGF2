#' plink QC 

PLINK=/home/g/gbader/spai/software/plink/plink
module load intel/15.0.6 R/3.2.3 

# --------------------------------------------------------
HM3=/scratch/g/gbader/spai/BaderLab/PNC/anno/hapmap3_pop/hg19/HM3_pops.hg19

# input genotype file
rawFile=/scratch/g/gbader/spai/NARSAD2014/input_files/CaseControlSNPchips/SUS19399_New/PLINK_140717_0448/SUS19399_New
outDir=/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399_New/hg19

### STOP. If you got the plink file from TCAG the coords don't match
### those in the SNP bpmap file from the manufacturer.
### First run fixPlink_coords.R in this folder.

outPref=${outDir}/SUS19399.hg19
# 0. convert to bed file
outDir_orig=/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399_New/hg38
mkdir -p $outDir_orig
outPref_orig=${outDir_orig}/SUS19399
$PLINK --map ${rawFile}.map --ped ${rawFile}.ped --make-bed --out $outPref_orig

# STOP here and run this script to convert coordinates to hg19
# before proceeding.
# ./plink_liftOver.sh

echo -e "\tPlot sample missingness" 
$PLINK --bfile ${outPref} --missing --out ${outPref}_imiss
### plot individual missingness in R
Rscript plink_plotMiss.R ${outPref}_imiss.imiss

echo -e "\tCompute outlying heterozygosity"
$PLINK --bfile ${outPref} --het --out ${outPref}_het
## plot in R
Rscript plink_plotHet.R ${outPref}_het.het

$PLINK --bfile $outPref --make-bed --out ${outPref}.sorted
## DON't UNCOMMENT THIS LINE
##outPref=${outPref}.sorted

echo "* Extract SNPs not in LD"
HIGHLD=high-LD-regions.txt
$PLINK --bfile ${outPref} --exclude $HIGHLD --indep-pairwise 50 5 0.2 --make-bed --out ${outPref}_prune

 echo "* Compute cryptic relatedness"
 Note prune.in contains the SNPs that were retained after LD filtering.
 prune.out contains those that were excluded by LD filtering.
time $PLINK --bfile ${outPref}_prune --extract ${outPref}_prune.prune.in --genome --out ${outPref}_IBS
Rscript plink_IBS.R ${outPref}_IBS.genome
echo "* Remove lines containing technical reps"
grep -v "\-R[123]" ${outDir}/fail-IBD.txt > ${outDir}/fail-IBD-clean.txt

# ---------------------------------------------------
# Section that combines with HapMap3 for ancestry ascertainment

echo "**** Computing divergent ancestry"
echo ""
time $PLINK --bfile ${outPref} --extract ${outPref}_prune.prune.in --make-bed --out ${outPref}_pruned
echo "  * Extract our SNPs from HM3 data"
$PLINK --bfile $HM3 --extract ${outPref}_prune.prune.in --filter-founders --make-bed --out ${outDir}/HapMapMini

echo -e "\t* Combine our data with HM3"
$PLINK --bfile ${outPref}_pruned --bmerge ${outDir}/HapMapMini.bed ${outDir}/HapMapMini.bim ${outDir}/HapMapMini.fam --geno 0.03 --make-bed --out ${outDir}/4Flips

# resolve unaligned SNPs by first getting a list of symmetric SNPs
echo "  * Compiling symmetric snps"
symSNP=${outDir}/syms.txt
ln -s ${outPref}_pruned.bim ${outDir}/top.bim
cat /dev/null > $symSNP;
awk '($5=="C" && $6 == "G"){print $2}' ${outDir}/top.bim >> $symSNP;
awk '($5=="G" && $6 == "G"){print $2}' ${outDir}/top.bim >> $symSNP;
awk '($5=="A" && $6 == "T"){print $2}' ${outDir}/top.bim >> $symSNP;
awk '($5=="T" && $6 == "A"){print $2}' ${outDir}/top.bim >> $symSNP;
echo "  * Merge using flip for missing snps"
$PLINK --bfile  ${outPref}_pruned --flip ${outDir}/4Flips-merge.missnp --make-bed --out ${outPref}_Align

echo "	* Merge with HapMapv3"
$PLINK --bfile ${outPref}_Align --bmerge ${outDir}/HapMapMini.bed ${outDir}/HapMapMini.bim ${outDir}/HapMapMini.fam --geno 0.03 --make-bed --out ${outDir}/HMDATA

echo "	* Generate pairwise IBS with HapMapv3"
$PLINK --bfile ${outDir}/HMDATA --genome --out ${outDir}/IBS_HM_DATA

$PLINK --bfile ${outDir}/HMDATA --read-genome ${outDir}/IBS_HM_DATA.genome --cluster --mds-plot 12 --out ${outDir}/MDS
###Rscript plink_plotMDS.R ${outPref}_pruned.fam ${outDir}/MDS.mds $outPref
######
# Exclude SNPs with >1% missing
# Exclude samples with 
cat ${outDir}/fail-het-qc.txt ${outDir}/fail-IBD-clean.txt | sort -k1 | uniq > ${outDir}/fail_inds.txt
echo "* Total failed samples"
wc -l ${outDir}/fail_inds.txt

echo "* Removing all markers and samples failing QC"
$PLINK --bfile ${outPref} \
    --maf 0.05 --geno 0.01 --hwe 0.000001 \
    --remove ${outDir}/fail_inds.txt --mind 0.1 --make-bed -out ${outPref}_-clean

echo "Run PCA to get ethnicity"
$PLINK --bfile ${outDir}/HMDATA --pca --out ${outDir}/HMDATA_PCA

###Extract just CEU-like
$PLINK --bfile ${outPref}_-clean --keep ${outDir}/ethnicities/CEU.txt --make-bed -out ${outPref}_-clean_CEU
$PLINK --bfile ${outPref}_-clean_CEU --pca --out ${outPref}_-clean_CEU_pca
