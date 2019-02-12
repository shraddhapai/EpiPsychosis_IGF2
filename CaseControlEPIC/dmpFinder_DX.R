#' find case/control DMPs
#' 
#' @details started with minfi::dmpFinder, but modified to include 
#' covariates in the null model
#' @param dat A MethylSet or a matrix.
#' @param outPref (char) output file prefix
#' @param excludeTechReps (logical) if TRUE, removes technical replicates
#' if FALSE, uses limma::duplicateCorrelation() to treat them as 
#' a blocking factor
#' @import limma
dmpFinder_DX <- function (dat, outPref,topX=0.5, excludeTechReps=FALSE,
						  incGeneticPCs=TRUE,incSlide=TRUE,incSex=TRUE,
						incPMI=TRUE)
{
	dt <- format(Sys.Date(),"%y%m%d")
	require(limma)
    M <- getBeta(dat)
	v	<- getVariance(M)
	idx <- which(v >= quantile(v,topX))
	cat(sprintf("\t Limiting to %i probes\n", length(idx)))

	dat <- dat[idx,]

	M	<- getBeta(dat)
	pd 	<- pData(dat)

	# build a single model with the covariates and contrast
	# of interest
	# recommended by Gordon Smyth as the 'acceptable method' in this
	# BioC support post
	# https://support.bioconductor.org/p/66251/
	# example >>>
	# fit <- eBayes(lmFit(data, model.matrix(~covariate+group)))
	# y <- topTable(fit, coef=3, number=Inf)
	# DEgenes <- y[y$adj.P.Val<0.05,]
	# <<<
	pd$AGE	<- as.integer(pd$AGE)
	pd$DX	<- factor(pd$DX,levels=c("control","case")) # control=1
	pd$SEX	<- factor(pd$SEX)
	pd$Slide <- factor(pd$Slide)
	pd$PMI <- as.numeric(pd$PMI)

	idx <- which(is.na(pd$PMI))
	if (any(idx)) {
		cat("Samples should all have PMI at this point!"); browser()
		cat(sprintf("Excluding the %i samples with na PMI\n", length(idx)))
		pd <- pd[-idx,]
		M	<- M[,-idx]
	}

	if (incGeneticPCs) {
		cat("* Including genetic PCs\n")
		if (incSlide) {
			cat("* Including microarray slide\n")
			if (!incSex) { 
				cat("this option hasn't been implementation\n") 
				browser()
			}
			design	<- model.matrix(~1+DX+SEX+AGE+Slide+C1+C2+PMI,data=pd)
		} else {
			cat("* Excluding microarray slide\n")
			if (incSex) {
				design	<- model.matrix(~1+DX+SEX+AGE+C1+C2+PMI,data=pd)
			} else {
				if (incPMI) {
					cat("* Excluding sex, but including PMI\n")
					design	<- model.matrix(~1+DX+AGE+C1+C2+PMI,data=pd)
				} else {
					cat("* Excluding sex and PMI\n")
					design	<- model.matrix(~1+DX+AGE+C1+C2,data=pd)
				}
			}
		}
	} else {
		cat("* EXCLUDING genetic PCs\n")
			if (!incSex) { 
				cat("this option hasn't been implementation\n") 
				browser()
			}
		if (incSlide) {
			cat("including microarray slide\n")
			design	<- model.matrix(~1+DX+SEX+AGE+Slide+PMI,data=pd)
		} else {
			cat("excluding microarray slide\n")
			design	<- model.matrix(~1+DX+SEX+AGE+PMI,data=pd)
		}
	}


	cat("Preview of model matrix\n")
	cat("-------------------------------------\n")
	print(head(design))
	cat("-------------------------------------\n")
	
	if (excludeTechReps) {
		idx <- grep("-R[23456]",pd$Sample_Name)
		cat(sprintf("* Excluding %i technical replicates\n",length(idx)))
		M <- M[,-idx]
		design <- design[-idx,]
		cat("\n\t* Calling lmFit\n")
    	fit <- eBayes(lmFit(M, design))
	} else {
		#  add technical replicate as a random factor
		cat("\n\t* Adding patient as a random factor\n")
		idx <- length(grep("-R[23456]", pd$Sample_Name))
		cat(sprintf("\t\t%i technical replicates identified\n",idx))
		sampName <- sub("-R[123456]","",pd$Sample_Name)

		cat(sprintf("%i unique individuals\n", length(unique(sampName))))
		tmp <- pd[!duplicated(pd$ID),]
		print(table(tmp[,c("DIST.DX","SEX")]))

		cat("\n\t* Calling lmFit: estimating correlation\n")
		corfit <- duplicateCorrelation(M,block=factor(sampName))
		save(corfit,file=sprintf("%s.corFit.Rdata",outPref))
		cat("Loading correlation\n")
		load(sprintf("%s.corFit.Rdata",outPref))
		print(table(sampName))

    	fit <- eBayes(lmFit(M, design,correlation=corfit$consensus))
	}
	
	cat("* Writing top tables and plotting nominal p\n")
	outFile <- sprintf("%s.dmp_nominalP.%s.pdf",outPref,dt)
	pdf(outFile)
	tryCatch({
	par(mfrow=c(2,2))
	out <- list()
	for (gp in c("DX")) {
		idx <- grep(gp, colnames(fit$coefficients))
		tmp <- colnames(fit$coef)[idx]
		cat(sprintf("\t%s\n",tmp))
		tt	<- topTable(fit,coef=tmp,number=Inf)


		hist(tt$P.Value,n=100,main=sprintf("%s: nominal p",tmp))
		tt <- tt[order(tt$adj.P.Val),]
		sig <- sum(tt$adj.P.Val<0.05)
		cat(sprintf("\t\tQ<0.05 = %i of %i (%i%%)\n", 
					sig,nrow(tt), round((sig/nrow(tt))*100)))

		curOut <- sprintf("%s.dmp_%s.topTable.%s.txt",outPref,tmp,dt)
		write.table(tt,file=curOut,sep="\t",
					col=TRUE,row=TRUE,quote=FALSE)
		out[[gp]] <- tt
	}
	},error=function(ex) {
		print(ex)
	},finally={
		dev.off()
	})
	
	out
}
