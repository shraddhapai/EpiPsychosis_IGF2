#' Runs blockFinder() with DX
#' 
#' @details started with minfi::dmpFinder, but modified to include 
#' covariates in the null model
#' @param dat A MethylSet or a matrix.
#' @param outPref (char) output file prefix
#' @param excludeTechReps (logical) if TRUE, removes technical replicates
#' if FALSE, uses limma::duplicateCorrelation() to treat them as a blocking
#' factor
blockFinder_DX <- function (dat, outPref,topX=0.5, excludeTechReps=TRUE) 
{
	dt <- format(Sys.Date(),"%y%m%d")
	require(limma)
    M <- getBeta(dat)
	v	<- getVariance(M)
	idx <- which(v >= quantile(v,topX))
	cat(sprintf("\t Limiting to %i probes\n", length(idx)))

	dat <- dat[idx,]
	pd	<- pData(dat)
	design <- model.matrix(~1+factor(DX)+factor(SEX)+factor(Slide),data=pd)
	design <- cbind(design, AGE=as.integer(pd$AGE))
	# coefs:2=DX; 3=SEX,4=Array,5=age
	
	###if (excludeTechReps) {
	idx <- grep("-R[23]",pd$Sample_Name)
	cat(sprintf("* Excluding %i technical replicates\n",length(idx)))
	design <- design[-idx,]
	dat <- dat[,-idx]
	print(dim(dat))

	browser()
	sink("blockFinder.log")
	tryCatch({
		cat(sprintf("Started at: %s", format(Sys.time())))
		warning("need to first make clusters")
		out <- blockFinder(MSet.genome,design,coef=2,what="Beta",
			nullMethod="bootstrap")
		browser()
	},error=function(ex) {
		print(ex)
	},finally={
		cat(sprintf("Ended at: %s", format(Sys.time())))
		sink(NULL)
	})
	save(out,file="blockFinder_output.Rdata")
###	} else {
###		# also add technical replicate as a random factor
###		cat("\n\t* Adding patient as a random factor\n")
###		cat(sprintf("\t\tHave at least %i technical replicates\n",
###			length(grep("-R1",pd$Sample_Name))))
###		sampName <- sub("-R[123]","",pd$Sample_Name)
###		cat("\n\t* Calling lmFit\n")
###		corfit <- duplicateCorrelation(M,design, block=factor(sampName))
###		save(corfit,file=sprintf("%s.corFit.Rdata",outPref))
###    	fit <- eBayes(lmFit(M, design,correlation=corfit$consensus))
######	}
###	
###	cat("* Writing top tables and plotting nominal p\n")
###	outFile <- sprintf("%s.dmp_nominalP.%s.pdf",outPref,dt)
###	pdf(outFile)
###	tryCatch({
###	par(mfrow=c(2,2))
###	out <- list()
###	for (gp in c("AGE","SEX","DX")) {
###		idx <- grep(gp, colnames(fit$coefficients))
###		tmp <- colnames(fit$coef)[idx]
###		cat(sprintf("\t%s\n",tmp))
###		tt	<- topTable(fit,coef=tmp,number=Inf)
###
###		hist(tt$P.Value,n=100,main=sprintf("%s: nominal p",tmp))
###		tt <- tt[order(tt$adj.P.Val),]
###		sig <- sum(tt$adj.P.Val<0.05)
###		cat(sprintf("\t\tQ<0.05 = %i of %i (%i%%)\n", 
###					sig,nrow(tt), round((sig/nrow(tt))*100)))
###
###		curOut <- sprintf("%s.dmp_%s.topTable.%s.txt",outPref,tmp,dt)
###		write.table(tt,file=curOut,sep="\t",col=TRUE,row=TRUE,quote=FALSE)
###		out[[gp]] <- tt
###	}
###	},error=function(ex) {
###		print(ex)
###	},finally={
###		dev.off()
###	})
###	
###	out
}
