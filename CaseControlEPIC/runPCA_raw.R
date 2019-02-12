#' run MDS on raw data.
require(minfi)
require(dataExplore) # getVar and pcaColored()

inFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/preprocessing/rawExp_161216.Rdata"
dt <- format(Sys.Date(),"%y%m%d")
outDir <- dirname(inFile)
nr <- ceiling(length(groupCols)/2)
outF <- sprintf("%s/PCA_raw_%s.pdf",outDir,dt)
topK <- 100*1000
#-------------
cat("* Loading raw data\n")
load(inFile)
cat("* Getting variance\n")
M 	<- getBeta(epicData)
v	<- getVariance(M)
v2	<- sort(v,decreasing=TRUE)
idx <- which(v >= v2[topK]) 
cat(sprintf("\t Limiting to %i probes\n", length(idx)))

pd <- pData(epicData)
colnames(M) <- pd$Sample_Name
M	<- M[idx,]

pdf(outF,width=11,height=11)
tryCatch({
	dataExplore::pcaColored(M,pd,groupToPlot=c("Cell.Group","SEX"),cex=1.5,
							las=1,bty='n',cex.axis=1.2,cex.lab=1.2) 
},error=function(ex){
	print(ex)
},finally={
	dev.off()
})

