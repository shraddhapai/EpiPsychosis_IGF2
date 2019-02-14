#' plot group means for identified DMPs
rm(list=ls())

rootDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398"
inFile <- sprintf("%s/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",
	rootDir)

outDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/agesexPMI_PC12_171129"

require(minfi)
cat("Loading data\n")
load(inFile)

probeList <- c("cg07096953","cg02613624","cg22956483","cg26401390")

# remove technical replicates (all but first)
pd <- pData(MSet.genome)
idx <- grep("-R[23456]",pd$Sample_Name)
MSet.genome <- MSet.genome[,-idx]; rm(pd)

# get order right in boxplot
pd <- pData(MSet.genome)
pd$DIST.DX <- factor(pd$DIST.DX,levels=c("Control","Schizophrenia","Bipolar"))
pData(MSet.genome) <- pd

print(table(pd$DIST.DX))

require(IlluminaEPICtools)
outFile <- sprintf("%s/IGF2probes.pdf",outDir)
pdf(outFile,width=8,height=8)
dmp_plotProbes_groups(MSet.genome,group="DIST.DX",baseline_group="control",
	group2="SEX",selProbes=probeList,plot_nr=2,plot_nc=2)
dev.off()

