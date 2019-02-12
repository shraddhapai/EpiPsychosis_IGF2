#' merge genetic-epigenetic
#'
#' @details Note that QC identified sample 90 being a technical 
# replicate of 95. This code also changes the genetic ID coding 
#' to reflect that.
#' @return mset, with first 5 genetic PCs added to phenotype table.
#' only samples with both genetic and epigenetic data will be retained
mergeWithGenetic <- function(mset,snpPCfile) {

	cat("* Adding genetic PCs\n")
	mds <- read.table(snpPCFile,as.is=T,h=F)
	
	colnames(mds)[1:ncol(mds)] <- c("FID","IID",
		paste("C",1:(ncol(mds)-2),sep=""))
	if (any(is.na(mds$IID))) mds <- mds[-grep("NA",mds$IID),] # exclude HapMap samples

	# recode sample 90 as 95 because these are technical replicates
	mds$IID[which(mds$IID == "90-R1")] <- "95-R4"
	mds$IID[which(mds$IID == "90-R2")] <- "95-R5"
	mds$IID[which(mds$IID == "90-R3")] <- "95-R6"

	mds <- cbind(mds, Sample_ID=sub("-R[123456]","",mds$IID))
	mds$Sample_ID <- as.character(mds$Sample_ID)
	
	pd <- pData(mset)
	midx <- match(pd$Sample_ID,mds$Sample_ID)
	if (all.equal(mds$Sample_ID[midx],pd$Sample_ID)!=TRUE) {
		cat("some methylation samples are missing genetic data\n")
		browser()

		cat("Exclude these\n")
		idx 	<- which(is.na(midx))
		mset	<- mset[,-idx]
		pd		<- pData(mset)  
		midx <- match(pd$Sample_ID,mds$Sample_ID)
		if (all.equal(mds$Sample_ID[midx],pd$Sample_ID)!=TRUE) {
			cat("even after excluding samples, genetic/methylation don't match!");
			browser()
		} else {
			cat("\tproblem resolved after excluding samples\n")
		}
	}
	mds <- mds[midx,]
	pd$C1 <- mds$C1
	pd$C2 <- mds$C2
	pd$C3 <- mds$C3
	pd$C4 <- mds$C4
	pd$C5 <- mds$C5
	pd$C6 <- mds$C6
	pData(mset) <- pd

	cat("After merging with genetic data, mset dimensions are:\n")
	print(dim(mset))

	mset
}
