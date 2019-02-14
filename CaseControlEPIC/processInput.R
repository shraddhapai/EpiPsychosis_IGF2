#' Preprocess EPIC pilot arrays
#'
#' instructions to adapt to methylationEPIC array are here:
#' https://bitbucket.org/hansenlab/illumina_epic
rm(list=ls())
require(minfi)

useQuantile <- TRUE ##### set to TRUE if quantile normalization is 
					##### desired

rootDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014"
baseIn <- sprintf("%s/input_files/CaseControlEPIC",rootDir)
baseOut <- sprintf("%s/output_files/CaseControlEPIC",rootDir)
outDir <- sprintf("%s/SUS19398/preprocessing",baseOut)
phenoFile <- sprintf("%s/NARSAD_sampleKey_160922_v2.txt",baseIn)

if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)
dt <- format(Sys.Date(),"%y%m%d")

# --------------------------------------------------------------
# Read raw data and perform QC
# --------------------------------------------------------------
logFile <- sprintf("%s/preprocess_readExp_QC_%s.log",outDir,dt)
sink(logFile,split=TRUE)
 tryCatch({
	# read pheno
	cat("\n* Reading pheno\n")
	sheet	<- read.metharray.sheet(sprintf("%s/SUS19398",baseIn))
	sampID	<- sub("-R[123]","",sheet$Sample_Name)
	sheet$Sample_ID <- sampID

	pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
	x <- merge(x=sheet,y=pheno,by="Sample_ID")
	if (nrow(x) < nrow(sheet)) {
		cat(sprintf("Only %i of %i samples matched, some unmatched. Check?\n",
		nrow(x),nrow(sheet)))
		browser()
	} else {
		cat(sprintf("%i of %i samples identified!\n", nrow(x),nrow(sheet)))
	}
	pheno <- x; rm(x)

	cat("* Clean pheno matrix\n")
	cat("------------------\n")
	cat("Add cell group\n")
	cat("------------------\n")
	tis <- rep("NeuN+", nrow(pheno))
	tis[which(pheno$Sample_Name%in% c("31","42"))] <- "NeuN-"
	pheno$Cell.Group <- tis
	print(table(pheno$Cell.Group))

	cat("Clean DIST.DX and add DX column\n")
	cat("------------------\n")
	print(table(pheno$DIST.DX))
	dx <- pheno$DIST.DX
	dx[dx %in% c("Schizophrenia","Bipolar")] <- "case"
	pheno$DX <- tolower(dx)
	print(table(pheno$DX))

	cat("Check other demographic variables\n-------------------\n")
	cat("Sex\n-----------\n")
	print(table(pheno$SEX))
	cat("Age\n-----------\n")
	print(summary(pheno$AGE))
	cat("Race\n----------\n")
	print(table(pheno$RACE))

	cat("Click enter if happy with pheno settings, else fix and rerun\n")
	browser()
				
	cat("\n* Reading experiment\n")
	t0 <- Sys.time()
	epicData<- read.metharray.exp(targets=pheno)
	t1 <- Sys.time(); print(t1-t0)

	cat("\tSaving raw data\n")
	save(epicData,file=sprintf("%s/rawExp_%s.Rdata",outDir,dt))

	
	cat("\n* Performing QC\n")
	pdfFile <- sprintf("%s/QualityControl_Raw_%s.pdf",outDir,dt)
	pd <- pData(phenoData(epicData))
	t0 <- Sys.time()
	qcReport(epicData, sampNames=pd$Sample_Name,
			 sampGroups=pd$age_group,pdf=pdfFile)
	t1 <- Sys.time()
	print(t1-t0)

	# exclude probes that fail detection filter
	cat("* Saving probes that failed detection\n")
	detP		<-detectionP(epicData)
	failed.01<- detP>0.01
	save(failed.01, file =sprintf("%s/failedProbes_%s.Rdata",outDir,dt))

	cat("\n* QC: Plot beta distribution of all samples\n")
	require(IlluminaEPICtools)
	pdf(sprintf("%s/betas_rawData_%s.pdf",outDir,dt))
	tryCatch({
		densityPlot_batch(epicData,8)
	},error=function(ex) {
		print(ex)
	},finally={
		dev.off()
	})
},error=function(ex){
	print(ex)
}, finally={
	cat("Closing log.\n")
	sink(NULL)
})

cat("\n* MDS plot\n")
	groupCols <- c("Array","Slide", "RACE", "SEX",
				   "Source","Cell.Group","DIST.DX")
	nr <- ceiling(length(groupCols)/2)
	
	cat("\n* QC: Multidimensional scaling by pheno groups\n")
	betas <- getBeta(epicData)
	pd <- pData(epicData)
	
	outF <- sprintf("%s/MDS_raw_%s.pdf",outDir,dt)
	pdf(outF,width=11,height=11)
	tryCatch({
		par(mfrow=c(1,1))
		k <- 100000
		for (g in groupCols) {
			mdsPlot(betas, numPositions=k,
				sampGroups=pd[,g],
				sampNames=pd$Sample_Name,
				main=sprintf("MDS on most variable %i pos\ncoloured by %s",
					k,g))
		}
	},error=function(ex){
		print(ex)
	},finally={
		dev.off()
	})


#------------------------------------------------------
# Phase 2: normalization
#------------------------------------------------------

logFile <- sprintf("%s/preprocess_NormalizeAndQC_%s.log",outDir,dt)
sink(logFile,split=TRUE)
rawFile <- sprintf("%s/rawExp_%s.Rdata",outDir,dt)
load(rawFile)

set.seed(42) # make reproducible

tryCatch({
	print(Sys.Date())
	print(Sys.time())
	if (useQuantile) {
		cat("Running noob\n")
		print(system.time(MSet.noob <- preprocessNoob(epicData)))

		cat("Quantile normalizing -- for this, excluding NeuN-\n")
		pd <- pData(MSet.noob)
		idx <- which(pd$Cell.Group %in% "NeuN-")
		MSet.noob <- MSet.noob[,-idx]
		print(table(pData(MSet.noob)$Cell.Group))
		# without this next line, noob output cannot be fed into quantile
		# https://github.com/kasperdanielhansen/minfi/issues/30
		MSet.noob@preprocessMethod <- c(rg.norm = 
										sprintf("Noob, dyeCorr=%s", TRUE))
		print(system.time(MSet.Qtl <- preprocessQuantile(MSet.noob)))
		MSet.genome <- mapToGenome(MSet.Qtl)
		save(MSet.Qtl,
		 	file=sprintf("%s/CaseControllEPIC_minfiOut_Qtl_%s.Rdata",
						  outDir,dt))
		save(MSet.genome,file=
			 sprintf("%s/CaseControlEPIC_MSetgenome_Qtl_%s.Rdata",
					 outDir,dt))
	} 

	cat("\n* MDS plot\n")
	groupCols <- c("Array","Slide", "RACE", "SEX",
				   "Source","Cell.Group","DIST.DX")
	nr <- ceiling(length(groupCols)/2)
	
	cat("\n* QC: Multidimensional scaling by pheno groups\n")
	betas <- getBeta(MSet.genome)
	pd <- pData(MSet.genome)
	
	outF <- sprintf("%s/MDS_normalized_%s.pdf",outDir,dt)
	pdf(outF,width=11,height=11)
	tryCatch({
		par(mfrow=c(1,1))
		k <- 100000
		for (g in groupCols) {
			mdsPlot(betas, numPositions=k,
				sampGroups=pd[,g],
				sampNames=pd$Sample_Name,
				main=sprintf("MDS on most variable %i pos\ncoloured by %s",
					k,g))
		}
	},error=function(ex){
		print(ex)
	},finally={
		dev.off()
	})

	sessionInfo()
},error=function(ex){ 
	print(ex)
},finally={
	cat("\n\nClosing log\n")
	print(Sys.time())
	sink(NULL)
})
