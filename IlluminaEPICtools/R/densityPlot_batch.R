#' density plot beta in chunks of samples
#'
#' @param rgset (RGChannelSet) raw experiment 
#' @param chunkSize (integer) num samples to plot on each page
#' @param verbose (logical) print reassuring messages
#' @importFrom minfi densityPlot
#' @import RColorBrewer
#' @return Nothing. Side effect of calling minfi::densityPlot() on each
#' subset of samples
#' @export
densityPlot_batch <- function(rgset, chunkSize=8,verbose=TRUE) {
	pal <- brewer.pal(n=chunkSize, name="Dark2")
	pd <- pData(rgset)
	
	numsamp <- ncol(rgset)
	maxpages <- ceiling(numsamp/chunkSize)
	ctr <- 1
	for (sidx in seq(1,numsamp,chunkSize)) {
		eidx <- min(sidx+(chunkSize-1),numsamp)
	if (verbose) cat(sprintf("\t%i-%i\n",sidx,eidx))
		densityPlot(rgset[,sidx:eidx],pd$Sample_Name[sidx:eidx],pal=pal)
		title(sprintf("Betas: Sample %i - %i\npage %i of %i",sidx,eidx,
					  ctr,maxpages))
		ctr <- ctr+1
	}
}
