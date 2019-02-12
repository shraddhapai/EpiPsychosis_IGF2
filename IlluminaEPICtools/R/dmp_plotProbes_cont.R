#' plot probe effects
#' 
#' @param mset (GenomicRatioSet)
#' @param xvar (char) which pheno to plot on x-axis
#' @param selIdx (integer) probe indices. these probes will be plotted
#' @param showSamp (integer) max number of probes to show
#' @param setSeed (integer) RNG seed
#' @param vals (char) column name in elementMetadata(mset) to plot
#' as part of title
#' @param ... params for plot() function
#' @export
dmp_plotProbes_cont <- function(mset,xvar="age",selIdx,showSamp=12L,
		setSeed=42L,vals,...){
	set.seed(setSeed)
	idx <- sample(selIdx, min(showSamp,length(selIdx)),F)

	betas <- getBeta(mset)
	pd <- pData(mset)[,xvar]
	loc <- getLocations(mset)
	vals <- elementMetadata(mset)[,vals]

	par(mfrow=c(4,3),bty='n',las=1)
	for (i in idx) {
		plot(pd,betas[i,],pch=16,cex=1.3,...)
		x <- summary(lm(betas[i,]~pd))
		abline(a=x$coef[1,1],b=x$coef[2,1],lwd=3,col='orange')
		title(sprintf("%s (%1.2e)",names(loc)[i],vals[i]))
	}
}
