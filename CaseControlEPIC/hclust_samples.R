#' unsupervised clustering of samples
#' 
#' @param mset
#' @param outFile - pdf for hclust
#' @param ... for IlluminaEpicTools::showDendro
hclustSamples <- function(mset,outFile,
	gps=c("Cell.Group","Slide","SEX","DX"),
	palName=list(Cell.Group="Dark2", Slide="Set3",SEX="Blues",DX="Dark2"),
	...) {
cat("* Preparing groups\n")

pd <- pData(mset)

for (g in gps) {
	pd[,g]		<- factor(pd[,g])
}

betas <- getBeta(mset)

# unsupervised clustering  
pdf(outFile,width=14,height=6)
tryCatch({
		cat("all samples\n")
IlluminaEPICtools::showDendro(M=betas,pheno=pd,groupPal=palName,
							  labRow=pd$Sample_Name,ttl="all",...)

},error=function(ex){
	print(ex)
},finally={
	dev.off()
})

}


