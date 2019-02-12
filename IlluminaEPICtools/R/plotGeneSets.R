#' plot trends for gene sets
#'
#' @param mset (MethylSet-Genomic)
#' @param gmt (char) path to gmt file
#' @param gene_GR (GRanges) 
plotGeneSets <- function(mset, gmt,gene_GR) {
	gs <- readPathways(gmt)
	p <- getLocations(mset)

	for (nm in gs) {
		g <- gene_GR[which(gene_GR$name2 %in% gs[[nm]])]
		cat(sprintf("%s: %i in list => %i genes\n",
			nm, length(gs[[nm]]), length(g)))
		plotAgeTrends(p,g)
	}
}
