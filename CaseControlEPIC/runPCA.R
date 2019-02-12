#' run PCA on beta values
require(dataExplore)
#'
#' @param mset
runPCA <- function(mset,outFile) {
	betas <- getBeta(mset)
	pd		<- pData(mset)
	pd$AGE	<- as.integer(pd$AGE)
	pd$PMI <- as.numeric(pd$PMI)
		
	.quantize <- function(x){
		qt <- rep("UNKNOWN", length(x))
		qt[which(x < quantile(x,0.33,na.rm=TRUE))] <- "LOW"
		qt[which(x >= quantile(x,0.33,na.rm=TRUE) & 
				 x < quantile(x,0.66,na.rm=TRUE))] <- "MED"
		qt[x>= quantile(x,0.66,na.rm=TRUE)] <- "HIGH"

		qt <- factor(qt,levels=c("LOW","MED","HIGH","UNKNOWN"))
	}

	pd$AGE_QT <- .quantize(pd$AGE)
	pd$PMI_QT <- .quantize(pd$PMI)

	pdf(outFile,width=9,height=9)
	tryCatch({
	pcaColored(betas,pd,
		groupToPlot=c("Cell.Group","Slide","DX","SEX",
					  "Source","AGE_QT","PMI_QT"))
	},error=function(ex) {
		print(ex)
	},finally={
		dev.off()
	})
}


