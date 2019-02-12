#' dendrogram from unsupervised clustering
#'
#' @param M (matrix) beta values
#' @param topVar (integer) top most-variable probes to show
#' @param pheno (data.frame) row order should match column order of M.
#' should contain column names for groupBy. Note that groupBy columns 
#' should be factors
#' @param groupPal (list) keys are columns in pheno by which data are to
#' be grouped. Values are RColorBrewer palette names
#' @param verbose (logical) print messages
#' @param ... parameters for squash::dendromat
#' @importFrom squash dendromat
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @return No value. Side effect of printing a dendrogram)
showDendro <- function(M, pheno, groupPal,topVar=10*1000L, verbose=TRUE,
					   ttl="",...) {

	if (verbose) cat("* preparing colours\n")
	clrmat <- matrix(NA,nrow=length(groupPal),ncol=nrow(pheno))
	rownames(clrmat) <- names(groupPal)
	ctr <- 1
	for (gpCol in names(groupPal)) {
		gpLev 	<- levels(pheno[,gpCol])
		pal		<- suppressWarnings(
			RColorBrewer::brewer.pal(
				name=groupPal[[gpCol]],n=max(3,length(gpLev)))
		)
	clrmat[ctr,]	<- pal[as.integer(pheno[,gpCol])]
	ctr <- ctr+1
	}
	clrmat <- clrmat[nrow(clrmat):1,]

	var 	<- IlluminaEPICtools::getVariance(M)
	M_top	<- M[order(var,decreasing=TRUE)[1:topVar],]
	d <- as.dist(1-cor(M_top))
	dendro	<- hclust(d,method="average")

	par(mar=c(6,6,4,3))
	squash::dendromat(dendro,t(clrmat),...)
	title(sprintf("%s: top %i\n",ttl,topVar))
}
