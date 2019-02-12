#' examine dynamics in regions of interest
#'
#' @details Finds subset of probes overlapping ROI. Then compares the rate
#' and direction of methylation change in this subset, relative to the 
#' backgound
#' @param mset (GenomicRatioSet) 
#' @param target_GR (GRanges) regions of interest 
#' @param xvar (char) pheno variable to use on x-axis
#' @param showSamp (integer) individual probes to show
#' @param ttl (char) string to put in plot title
#' @param ... param for sample_nonParametric()
#' @export
regionDynamics <- function(mset,target_GR,xvar="age",showSamp,ttl="",
	...) {
	loc <- getLocations(mset)
	betas <- getBeta(mset)
	pd <- pData(mset)
	cat(sprintf("\t%i target regions\n",length(target_GR)))

	require(limma)
	newdf <- data.frame(x=pd[,xvar])
	mod <- model.matrix(~1+x,data=newdf)
	fit <- lmFit(betas,mod)
	coef <- fit$coefficients
	idx <- findOverlaps(loc,target_GR);
	idx <- unique(idx@from)
	cat(sprintf("\t%i probes; %i samples\n",length(idx), ncol(betas)))

	# are they mostly increasing ? decreasing ? staying the same? 
	tmp <- list(in.region=coef[idx,2],overall=coef[,2])
	numReps <- 1e4
	pval <- sample_nonParametric(tmp[[1]],tmp[[2]],type="twotailed",
				R=numReps)
	cat(sprintf("\tForeground median = %1.2e (two tailed p = %1.2e)\n", 
				median(tmp[[1]]),pval))
	cat(sprintf("\tNum: Fore=%i ; Back = %i\n", length(tmp[[1]]),
				length(tmp[[2]])))

	par(las=1,bty='n',cex.axis=1.3)
	boxplot(tmp,pars=list(boxwex=0.4),
			main=sprintf("%s:%i regions (F=%i,B=%i)\nWMW p=%1.2e",
						 ttl,length(target_GR),
						 length(tmp[[1]]),length(tmp[[2]]),pval
						 ),
			ylab="Beta slope with age")
	abline(h=0,lwd=3,col='red')

	# plot most variable probes in this set
	v <- getVariance(betas[idx,]) # get 
	dmp_plotProbes(mset[idx,],xvar=xvar,selIdx=which(v > quantile(v,0.9)),
				   vals="qval")
	# show an average plot with mean and CI. are they hypermethylated/
	# hypomethylated? 
	
	# show distribution of variance relative to the background
	# are they unusually active compared to those not in this region?
	cat("\n")
	out <- list(num_regions=length(target_GR),
				size_fg=length(tmp[[1]]), size_bg=length(tmp[[2]]),
				slope_fg=median(tmp[[1]]),slope_fg=median(tmp[[2]]),
				pvalue=pval,numReps=numReps)
	out
}

