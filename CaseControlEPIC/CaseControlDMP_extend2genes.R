#' write DMP regions + enhancers/promoters in nearby genes, to
#' filter SNPs in region

rm(list=ls())
require(GenomicRanges)

# dmp hits to start from
inFile<- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/age_sex_c1c2_slide/dmp.blockTechReps.dropSNPs.dmp_factor(DX)control.topTable.161205.FDR0.30.locations.txt"
# only look at regions with DMPs that have Q < this threshold
FDRthresh <- 0.2

annoFiles <- list(
	enh1="/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/enhancer_defs/Yu_DLPFC_AllBrain_enhprom.bed",
	enh2="/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/enhancer_defs/Dennis_merged.bed",
	ctcf="/Users/shraddhapai/Google Drive/genome_annotation/ENCODE/GSM733765.bed",
	ctcf2="/Users/shraddhapai/Google Drive/genome_annotation/ENCODE/GSM1022662.bed")

dat <- read.delim(inFile,sep="\t",h=T,as.is=T)
dat <- subset(dat,qval < FDRthresh)
dat_GR <- GRanges(dat$seqnames,IRanges(dat$start,dat$end))
dat_GR$name <- rownames(dat)

# scan broader region
region_GR <- resize(dat_GR,width=20*1000,fix="center")
icr_GR <- GRanges("chr11",IRanges(1999227, 2218058))
icr_GR$name<- "H19-ICR"
region_GR <- c(region_GR,icr_GR)

cat("* Finding subset of gene_GR overlapping enhancers/CTCF\n")
reg_GR <- GRanges()
for (nm in names(annoFiles)) {
	cat(sprintf("%s\n",nm))
	
	dat <- read.delim(annoFiles[[nm]],sep="\t",h=F,as.is=T)
	gr <- GRanges(dat[,1],IRanges(dat[,2],dat[,3]))
	seqlevels(gr, force=TRUE) <- seqlevels(region_GR)

	x <- intersect(region_GR, gr)
	reg_GR <- c(reg_GR, x)
}
reg_GR <- reduce(reg_GR)
idx <- which(width(reg_GR)<1000)
if (any(idx)) {
	cat(sprintf("Resizing %i wins < 1000bp to 1000bp\n", length(idx)))
	reg_GR[idx] <- resize(reg_GR[idx],fix="center",width=1000)
}
reg_GR$name <- paste("CaseControl_Reg_",1:length(reg_GR),sep="")
cat(sprintf("%i ranges with width distribution:\n", length(reg_GR)))
print(summary(width(reg_GR)))

# keep:
# 1. regions in dat_GR -- DMP 2k window, reduced
# 2. segments of gene_GR that overlap regulatory elements
dat_GR <- resize(dat_GR,width=2000,fix="center")
dat_GR <- reduce(dat_GR)
dat_GR$name <- paste("CaseControl_DMP_",1:length(dat_GR),sep="")
keep_GR <- c(dat_GR, reg_GR)
keep_df <- as.data.frame(keep_GR)

options(scipen=10)
outDir	<- dirname(inFile)
dt <- format(Sys.Date(),"%y%m%d")
outFile	<- sprintf("%s/dmpCaseControl_Genes_%s.txt",outDir,dt)
cat(sprintf("Writing %i regions\n",nrow(keep_df)))
write.table(keep_df,file=outFile,sep="\t",col=FALSE,row=FALSE,quote=FALSE)




