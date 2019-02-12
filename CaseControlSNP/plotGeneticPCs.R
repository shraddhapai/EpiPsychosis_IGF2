#' plot case/control population against genetic pcs.
rm(list=ls())

rootDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files"
cleanDat <- sprintf("%s/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171029.Rdata",
	rootDir)
outDir <- sprintf("%s/CaseControlSNP/plinkQC",rootDir)

load(cleanDat)

cat("\n* Preparing for DMP analysis\n")
pd <- as.data.frame(pData(MSet.genome))
pd <- pd[!duplicated(pd$ID),]
pd <- pd[,c("DX","C1","C2","C3","C4","C5","C6")]

require(GGally)
pdf(sprintf("%s/caseControl_geneticPCs.pdf",outDir),width=11,height=11)
ggpairs(pd,aes(colour=DX,alpha=0.4))
dev.off()
