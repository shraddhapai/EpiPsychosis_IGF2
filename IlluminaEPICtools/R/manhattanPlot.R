#' Manhattan plot of probe p results
#'
#' @param mset (MSet.genome) 
#' @param dmpTable (data.frame) dmp results. Output of topTable()
#' function, must contain probe IDs as rownames, and the columns
#' "pval" (nominal p value) and "qval" (FDR corrected Q)
#' @export
manhattanPlot <- function(mset, dmpTable) {
	locs <- as.data.frame(getLocations(mset))
	locs$probe <- rownames(locs)
	dmp$probe <- rownames(dmp)
	x <- merge(x=locs,y=dmp,by="probe")
	x <- x[,c("start","seqnames","pval","probe")]
	x$seqnames <- sub("chr","",x$seqnames)

	idx <- which(x$seqnames %in% "X")
	if (any(idx)) x$seqnames[idx] <- 23
	idx <- which(x$seqnames %in% "Y")
	if (any(idx)) x$seqnames[idx] <- 24
	idx <- which(x$seqnames %in% "M")
	if (any(idx)) x$seqnames[idx] <- 25
	x$seqnames <- as.integer(x$seqnames)
	colnames(x) <- c("BP","CHR","P","SNP")

	highlight <- NULL
	idx <- which(dmp$qval < 0.05)
	if (any(idx)) highlight <- dmp$probe[idx]
 
	manhattan(x, 
		suggestiveline = F, 
		genomewideline = -log10(0.05/nrow(dmp)), 
		highlight=highlight,cex=0.7)
	#brewer.pal(5, "Set2") , cex=0.7)

}
