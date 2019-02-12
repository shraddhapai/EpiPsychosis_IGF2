#' for DMPs, gets mean delta beta for the groups for each hit.
rm(list=ls())
suppressMessages(suppressWarnings(require(minfi)))

rootDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398"
cleanDat <- sprintf("%s/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",rootDir)

dmpFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/agesexPMI_PC12_171129/dmp..dmp_DXcase.topTable.171129.txt"

cat("* reading dmps\n")
dmp <- read.delim(dmpFile,sep="\t",h=T,as.is=T)
cat("* loading methylation data\n")
load(cleanDat)

# remove tech reps
samps <- pData(MSet.genome)$Sample_Name
idx <- grep("-R[23456]",samps)
noTechReps <- MSet.genome[,-idx] # remove tech duplicates
pd <- pData(noTechReps)
tmp <- pd$DIST.DX
tmp[grep("Control",tmp)] <- "CONT"
tmp[grep("Schiz",tmp)] <- "SCZ"
tmp[grep("Bipolar",tmp)] <- "BPD"
pd$DX2 <- factor(tmp,levels=c("CONT","SCZ","BPD"))

probeNames <- c("cg07096953","cg02613624","cg22956483")

M <- getBeta(noTechReps)
igf2 <- M[probeNames,]

require(reshape2)
igf2 <- melt(igf2)

g <- unique(pd$DX2)
groups <- sapply(g,function(x) { rownames(pd)[pd$DX2==x]})
names(groups) <- g

for (k in probeNames){ 
	cur <- igf2[which(igf2$Var1 %in% k),]
	ct <- mean(cur$value[which(cur$Var2 %in% groups[["CONT"]])])
	sz <- mean(cur$value[which(cur$Var2 %in% groups[["SCZ"]])])
	bd <- mean(cur$value[which(cur$Var2 %in% groups[["BPD"]])])
	both <- mean(cur$value[which(cur$Var2 %in% c(groups[["BPD"]],groups[["SCZ"]]))])
	
	cat(sprintf("%s: CONT = %1.2f\n", k, ct))
	cat(sprintf("\tCONT - SCZ (%1.2f) = %1.2f\n", sz, ct-sz))
	cat(sprintf("\tCONT - BPD (%1.2f) = %1.2f\n", bd, ct-bd))
	cat(sprintf("\tCONT - BOTH (%1.2f) = %1.2f\n", both, ct-both))
	cat("\n")
}
