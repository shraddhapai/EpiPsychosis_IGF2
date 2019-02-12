#' plot age effect for a set of genomic intervals
#'
#' @param mset (MethylSet-genomic)
#' @param subject (char or GRanges) see getProbeDistr
#' @import gplots
#' @export
plotAgeTrends <- function(mset,subject) {
	locs	<- getLocations(mset)

	cat("* Average betas by age\n")
	p	<- pData(mset)
	gp	<- factor(p$age_group,levels=1:5)
	betas	<- getBeta(mset)
	b2 <- matrix(NA, nrow=nrow(betas),ncol=length(levels(gp)))
	colnames(b2) <- paste("age_", as.character(levels(gp)),sep="")
	ctr <- 1
	for (g in levels(gp)) {
		b2[,ctr] <- rowMeans(betas[,which(gp %in% g),drop=FALSE])
		cat(sprintf("group %s: %i samples\n", g, sum(gp %in% g)))
		ctr <- ctr+1
	}
	rownames(b2) <- rownames(betas)
	b2 <- b2/b2[,1]
	betas <- b2; rm(b2)

	cat("* Get probe overlap\n")
	ol 		<- getProbeDistr(locs, subject,getSubset=TRUE)

	dat <- ol[["overlap"]]
	for (nm in names(dat)) {
		cur	<- dat[[nm]]
		b <- betas[rownames(betas)%in%names(cur),]
		browser()
		
	}
}
