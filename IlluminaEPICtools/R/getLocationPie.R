#' pie chart of location breakdown
#' 
#' @details subtractive. Starts with first category and moves down from
#' there. Does not guarantee that categories are mutually exclusive.
#' @param gr (GRanges) ranges of interest
#' @param bedList (data.frame) names of categories (Type) and path to bed 
#' files with ranges (File)
#' @param verbose (logical) print messages
#' @import RColorBrewer
getLocationPie <- function(gr, bedList,pal="Spectral",verbose=TRUE) {
	# additional spot for residual
	props <- numeric(nrow(bedList)+1)
	names(props) <- c(bedList$Type,"other")
	# loop over bed list.
	for (k in 1:nrow(bedList)) {
		bed <- read.delim(bedList$File[k],sep="\t",as.is=T)
		bed_GR <- GRanges(bed[,1],IRanges(bed[,2],bed[,3]))
		seqlevels(bed_GR,force=TRUE) <- seqlevels(gr)

		ol <- findOverlaps(gr,bed_GR)
		ol <- unique(ol@from) # query hits
		props[k] <- length(ol)
		cat(sprintf("\t%s: %i overlaps\n", bedList$Type[k],length(ol)))
		gr <- gr[-ol] # remove from future counts
		k <- k+1
	}
	props[length(props)] <- length(gr)
	pct <- round((props/sum(props))*100)

	names(props) <- sprintf("%s (%s%%)", names(props),pct)

	# plot pie chart
	pie(props, clockwise=TRUE,
		col=RColorBrewer::brewer.pal(n=length(props),name=pal),
		cex=1)
}
