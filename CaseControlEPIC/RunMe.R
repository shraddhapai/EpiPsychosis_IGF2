# run case-control analyses
rm(list=ls())
require(IlluminaEPICtools)
require(minfi)
source("correctBatch.R")

dropSNPs	<- TRUE

# EPIC array data
rootDir 	<- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398"
inData		<- sprintf("%s/preprocessing/CaseControlEPIC_MSetgenome_Qtl_161216.Rdata",rootDir)
inDir		<- dirname(inData)
outRoot <- rootDir
failedProbes<- sprintf("%s/failedProbes_161216.Rdata",inDir)
xReactiveFile	<- "/Users/shraddhapai/Google Drive/genome_annotation/chip_annotation/Illumina_MethylationEPIC/McCartneyEvans_2016_GenomicsData_TableS2.txt"

# Genetic psych array data
snpPCFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlSNP/plinkQC/HMDATA_PCA.eigenvec"

annoFiles 	<- read.delim("Annotations.txt",sep="\t",h=T,as.is=T)

dt <- format(Sys.Date(),"%y%m%d")
cleanDat <- sprintf("%s/CaseControlEPIC_CLEAN_171029.Rdata",
					dirname(inData))
# ------------------------------------------------------------
# Cleaning the data
logFile <- sprintf("%s/cleaningData_%s.log",outRoot,dt)
sink(logFile,split=TRUE)
tryCatch({
	cat("* Loading data\n")
	print(system.time(load(inData)))
	
	cat("\n* Checking for mismatched sex\n")
	source("checkSex.R")
	flagSamps <- checkSex(MSet.genome)
	
	cat("* Dropping unreliable probes\n")
	# overlap polymorphisms
	if (dropSNPs) {
		cat("\tOverlap polymorphisms:")
		before <- dim(MSet.genome)
		MSet.genome <- dropLociWithSnps(MSet.genome,snps=c("CpG","SBE"),
							maf=0.05)
		after <- dim(MSet.genome)
		cat(sprintf(" %i probes\n",before[1]-after[1]))
	}
	
	# remove probes that fail detection
	load(failedProbes)
	failedProbes <- rownames(failed.01)[rowMeans(failed.01)>0.2]
	cat(sprintf("\tFailed detection: %i probes fail in >20%% samples\n",
		length(failedProbes)))
	badProbes <- failedProbes
	
	xReac <- read.delim(annoFiles[which(annoFiles[,1]%in%"EPIC_xReactive"),3],
			sep="\t",h=F,as.is=T)[,1]
	cat(sprintf("\tCross-reactive probes: %i \n",length(xReac)))
	badProbes <- c(badProbes,xReac)
	
	MSet.genome <- MSet.genome[-which(rownames(MSet.genome)%in%badProbes),]
	cat(sprintf("\tAfter removing all bad probes: %i probes left\n",
				dim(MSet.genome)[1]))
	
	cat("* Run PCA to identify samples to exclude\n")
	source("runPCA.R")
	outFile <- sprintf("%s/PCA_beforeBatchCorrect_%s.pdf",inDir,dt)
	runPCA(MSet.genome, outFile)
	
	cat("\n* hclust of samples\n")
	source("hclust_samples.R")
	outFile <- sprintf("%s/hClust_funnorm_%s.pdf",inDir,dt)
		hclustSamples(MSet.genome,outFile,topVar=100*1000L)

	cat("---------------------------\n")
	cat("QC\n")
	cat("---------------------------\n")
	samp2rm <- c("200590500059_R03C01")
	pd <- pData(MSet.genome)
	idx <- which(rownames(pd)%in% samp2rm)
	cat(sprintf("%i samples being excluded following PCA viz\n",
			length(idx)))
	idx2 <- which(pd$Sample_Name %in% c(78))
	cat(sprintf(
	"\t%i samples removed because genetic analysis flagged as duplicates\n",
			length(idx2)))
	idx <- unique(c(idx,idx2))
	cat(sprintf("Samples excluded: { %s }\n", 
		paste(pd$Sample_Name[idx],collapse=",")))
	MSet.genome <- MSet.genome[,-idx]; rm(pd)

	cat("\t* Recoding sample 90 as sample 95\n")
	pd <- pData(MSet.genome)
	idx <- grep("90",pd$Sample_ID)
	pd$Sample_ID[idx] <- "95"
	pd$Sample_Name[grep("90-R1", pd$Sample_Name)] <- "95-R4"
	pd$Sample_Name[grep("90-R2", pd$Sample_Name)] <- "95-R5"
	if (any(c(grep("90",pd$Sample_Name),grep("90",pd$Sample_ID)))) {
		cat("Sample 90 wasn't fully removed\n")
	}
	pData(MSet.genome) <- pd
	cat("Updated rows:\n")
	print(pd[idx,1:6])
	rm(pd)

	x <- dim(MSet.genome)
	cat(sprintf("After all QC, %i markers and %i samples left\n",
				x[1],x[2]))

	cat("* Merge with genetic data\n")
	source("mergeWithGenetic.R")
	MSet.genome <- mergeWithGenetic(MSet.genome,snpPCfile)

	# SVA and ComBaT batch correction
	cat("* Running SVA\n")
	#x <- correctBatch(MSet.genome,runCombat=FALSE)
	# Peter Langfelder recommends not to run ComBaT to correct batch
	# if running lmFit:
	# see https://support.bioconductor.org/p/80685/
	# https://support.bioconductor.org/p/69328/#69333

	cat("* Saving clean data for downstream analyses\n")
	save(MSet.genome,file=cleanDat)
},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})


