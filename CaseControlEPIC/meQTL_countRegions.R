#' take significant cis-meQTL results and identify how many regions 
#' of the 18 significant DMRs had this.
rm(list=ls())
require(minfi)

rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files"
hitFile <- sprintf("%s/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/meQTL/180110/regular_lmer_genoDisease0/withC1C2/eQTL.results.180110.txt",rootDir)
dmpFile <- sprintf("%s/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/agesexPMI_PC12_171129_genes.txt",
	   rootDir)
mData <- sprintf("%s/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",
				 rootDir)
cat("reading qtls\n")
qtls <- read.delim(hitFile,sep="\t",h=T,as.is=T)
cat("reading dmps\n")
dmps <- read.delim(dmpFile,sep="\t",h=T,as.is=T)
cat("loading mData\n")
load(mData)

dmps <- subset(dmps,z_sidak_p < 0.05)
qtls <- subset(qtls, FDR<0.05)
dmp_GR <- GRanges(dmps[,1],IRanges(dmps[,2],dmps[,3]))
locs <- getLocations(MSet.genome)
locs <- locs[unique(qtls$cpg)]

qtl_regions <- subsetByOverlaps(dmp_GR,locs)

cat(sprintf("%i DMP regions (%i CpGs) have significant cis-meQTLs; these involve %i unique SNPs\n",
	length(qtl_regions),length(locs), length(unique(qtls$snps))))




