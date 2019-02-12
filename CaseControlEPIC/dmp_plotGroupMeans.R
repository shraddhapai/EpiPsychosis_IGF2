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

####def <- par("mar")
####par(mfrow=c(4,6),mar=c(1,1,1,1))
###gp <- matrix(NA, nrow=nrow(betas),ncol=2)
###for (k in 1:nrow(betas)) {
###	x1 <- betas[k,which(pd$DX == "case")]
###	x2 <- betas[k,which(pd$DX == "control")]
###	gp[k,1] <- mean(x1)
###	gp[k,2] <- mean(x2)
####	boxplot(list(case=x1,control=x2),title=rownames(betas)[k],
####		ylim=c(0,1))
####	abline(h=0.5,lty=3,col='red')
###}
###
###
###par(mfrow=c(2,2)) #,mar=c(5,5,3,3))
###plot(gp[,1],gp[,2],xlab="case",ylab="control",
###main="mean group %M for top DMPs",col=rgb(0,0,0,0.4))
###abline(0,1,col='red',lty=3,lwd=2)
###hist(gp[,2], main="mean control % for DMR")
###df <- gp[,2]-gp[,1]
###hist(df,main="diff in group %M")


