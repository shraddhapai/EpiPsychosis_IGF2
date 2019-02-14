#' for DMPs, gets mean delta beta for the groups for each hit.
#' CEU males only.

rm(list=ls())
suppressMessages(suppressWarnings(require(minfi)))
require(ggplot2)
options(warn=4)

rootDir <- "~/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398"
cleanDat <- sprintf("%s/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",rootDir)
snpPCFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/SUS19399.hg19.sorted_-clean_CEU_pca.eigenvec"
updatedSampKey <- "/home/shraddhapai/Epigenetics/NARSAD/input_files/NARSAD_sampleKey_171129.txt"
dt <- format(Sys.Date(),"%y%m%d")

logFile <- sprintf("IGF2_probeDiff_CEUmale_%s.log",dt)
pdf(sprintf("IGF2_probeDiff_CEUmale_%s.pdf",dt))
sink(logFile,split=TRUE)

tryCatch({
	cat("* loading methylation data\n")
	load(cleanDat)
	pd <- pData(MSet.genome)
	
	# update PMI for last three samples
	p2 <- read.delim(updatedSampKey,sep="\t",h=T,as.is=T)
	pd$PMI <- suppressWarnings(as.numeric(pd$PMI))
	idx <- which(is.na(pd$PMI)) # find missing pmi
	cat("Filling in missing PMIs\n")
	for (k in idx) {
		idx2 <- which(as.character(p2$Internal.ID) == pd$Sample_ID[k])
		cat(sprintf("\tSample %s\n", pd$Sample_ID[k]))
		pd$PMI[k] <- p2$PMI[idx2]
	}
	if (TRUE) {
		cat("*********\n")
		cat("**** Limiting to CEU samples ****\n")
	
		pd <- pd[,-which(colnames(pd) %in% paste("C",1:6,sep=""))]
		pData(MSet.genome) <- pd
		source("mergeWithGenetic.R")
		MSet.genome <- mergeWithGenetic(MSet.genome,snpPCFile)
		cat("Merge with new successful!\n")
		pd <- pData(MSet.genome)
		idx <- which(pd$SEX %in% "Male")
		MSet.genome <- MSet.genome[,idx]
		pd <- pData(MSet.genome)
		cat(sprintf("*** FILTERING SEX=Male: %i samples left\n", nrow(pd)))
	}
	
	# remove tech reps
	samps <- pData(MSet.genome)$Sample_Name
	idx <- grep("-R[23456]",samps)
	noTechReps <- MSet.genome[,-idx] # remove tech duplicates
	pd <- pData(noTechReps)
	
	probeNames <- c("cg07096953","cg02613624","cg22956483","cg26401390")
	
	M <- getBeta(noTechReps)
	igf2 <- M[probeNames,]
	
	require(reshape2)
	igf2 <- melt(igf2)
	
	pd <- pData(noTechReps)
	pd <- pd[,c("DIST.DX","SEX","PMI","AGE","C1","C2")]
	pd$Var2 <- rownames(pd)
	pd$DX <- rep("Control",nrow(pd))
	pd$DX[which(pd$DIST.DX %in% c("Bipolar","Schizophrenia"))] <- "Case"
	
	igf2 <- merge(x=igf2,y=pd,by="Var2")
	
	for (k in probeNames){ 
		cur <- igf2[which(igf2$Var1 %in% k),]
	
		print(table(cur$DX))
		m1 <- lm(value~AGE+PMI+C1+C2,data=cur)
		m2 <- lm(value~DX+AGE+PMI+C1+C2,data=cur)
		res <- anova(m1,m2)
		pval <- res[["Pr(>F)"]]
		cat(sprintf("%s: pvalue from nested AOV = %1.2e\n",k,pval[2]))

	
		cur <- as.data.frame(cur)
		cur$value <- cur$value*100
		case_val <- cur$value[which(cur$DX %in% "Case")]
		ctrl_val <- cur$value[which(cur$DX %in% "Control")]
		cat(sprintf("%s: Case - Control = %1.2f\n", k, mean(case_val)-mean(ctrl_val)))
		p <- ggplot(cur,aes(x=DX,y=value)) + geom_boxplot(aes(fill=DX))
		p <- p + ggtitle(sprintf("%s (p= %1.2e)",k, pval[2]))
		p <- p + theme_bw() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
		print(p)
		cat("--------------------\n")
		cat("\n")
	}
}, error=function(ex) {
	print(ex)
}, finally={
	dev.off()
	sink(NULL)
})
