#' write GMT from results of pathway ORA
#' 
#' @details To be run when you have pathway-level p and Q values from
#' FET and binomial test, and you just want to generate the gmt to
#' visualize an EnrichmentMap
#' @param pathwayRes (data.frame) table output of dmp_pathwayORA. 
#' columns are: Pathway, n_fg, n_bg, Hypergeom_p, Binomial_p,
#' Hypergeom_Q, Binom_Q, genes_fg
#' @param pathwayList (list) keys are pathways, values are member genes
#' @param fg_GR (GRanges) foreground DMP probes. These will be used
#' to filter for genes that contribute to signal
#' @param gene_GR (GRanges) domains for each gene, may include upstream or
#' downstream of TSS.
#' @param geneDomain (integer) num kilobases to extend gene upstream of TSS
#' @param outFile (char) output file
#' @export
dmp_pathwayORA_writeGMTFromRes <- function(pathwayFile, selPathways, fg_GR,
	gene_GR,geneDomain=5,outFile){
	pathwayList <- readPathways(pathwayFile)
	# define gene domain
	gene_GR <- resize(gene_GR,width=width(gene_GR)+(geneDomain*1000),
	  fix="end")

	if (file.exists(outF)) unlink(outF)
	system(sprintf("touch %s",outF))
	for (nm in selPathways) {
		ol <- findOverlaps(fg_GR,gene_GR)
		curg <- unique(gene_GR$name[ol@to])
		g <- intersect(pathwayList[[nm]], curg)
		cat(sprintf("%s\t%s\t%s\n",nm,nm,
				paste(pathwayList[[nm]],collapse="\t")),file=outF,
				append=TRUE)
	}
}
