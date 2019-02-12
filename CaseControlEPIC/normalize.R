#------------------------------------------------------
# Phase 2: normalization
#------------------------------------------------------
require(minfi)

logFile <- sprintf("%s/preprocess_NormalizeAndQC_%s.log",outDir,dt)
sink(logFile,split=TRUE)
rawFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/preprocessing/rawExp_160922.Rdata"
load(rawFile)

set.seed(42) # make reproducible

tryCatch({
	print(Sys.Date())
	print(Sys.time())
	if (useQuantile) {
		cat("Running noob\n")
		print(system.time(MSet.noob <- preprocessNoob(epicData)))
		cat("* Fixing pheno\n")

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
	rm(pd)
	
	cat("\n* MDS plot\n")
	groupCols <- c("Array","Slide", "RACE", "SEX",
				   "Source","Cell.Group","DIST.DX")
	nr <- ceiling(length(groupCols)/2)
	
	cat("\n* QC: Multidimensional scaling by pheno groups\n")
	betas <- getBeta(MSet.genome)
	
	outF <- sprintf("%s/MDS_normalized_%s.pdf",outDir,dt)
	pd <- pData(MSet.genome)
	pdf(outF,width=11,height=11)
	tryCatch({
		par(mfrow=c(1,1))
		k <- 10000
		for (g in groupCols) {
			cat(sprintf("Group: %s\n",g))
			mdsPlot(betas, numPositions=k,
				sampGroups=pd[,g],
				sampNames=pd$Sample_Name,
				#xlim=c(-20,20),ylim=c(-10,10),
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
