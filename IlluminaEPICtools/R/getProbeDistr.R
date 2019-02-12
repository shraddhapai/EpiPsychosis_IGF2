#' get probe subset overlapping ROI. NOT strand-aware
#'
#' @param query (GRanges) probe locations
#' @param subject (char or GRanges) if char, vector of bed files to read in as 
#' GRanges and compute ol.
#' @param getSubset (logical) if TRUE, returns probes that overlap gene ranges
#' @return probe breakdown
#' @export
#' @import GenomicRanges
getProbeDistr <- function(query,subject,getSubset=FALSE) {
	cat(sprintf("Got %i bed files\n", length(subject)))
	out 	<- list()
	stats	<- numeric()
	if (class(subject)=="character") {
	for (f in subject) {
		baseF <- basename(f)
		cat(sprintf("%s\n",baseF))
		dat 	<- read.delim(f,h=F,as.is=T,sep="\t")
		dat_GR	<- GRanges(dat[,1],IRanges(dat[,2],dat[,3]))
		dat_GR$name <- dat[,4]
		seqlevels(dat_GR,force=TRUE) <- seqlevels(query)

		ol <- findOverlaps(query,dat_GR)
		ol	<- cbind(ol@queryHits,ol@subjectHits)
		if (getSubset) {
			tmp <- query[ol[,1]]
			tmp$REGION_START	<- start(dat_GR[ol[,2]])
			tmp$REGION_END		<- end(dat_GR[ol[,2]])
			tmp$REGION_NAME		<- dat_GR$name[ol[,2]]

			out[[baseF]] <- tmp
		}

		stats <- c(stats, length(unique(ol[,1])))
	}
	names(stats) <- basename(subject)

	} else {
		ol <- findOverlaps(query,subject)
		ol	<- cbind(ol@queryHits,ol@subjectHits)
		if (getSubset) {
			tmp <- query[ol[,1]]
			tmp$REGION_START	<- start(subject[ol[,2]])
			tmp$REGION_END		<- end(subject[ol[,2]])
			tmp$REGION_NAME		<- subject$name[ol[,2]]

			out[[1]] <- tmp
		}
		stats <- c(stats, length(unique(ol[,1])))
	}
	if (getSubset) {
		out <- list(stats=stats,overlap=out)
	} else {
		out <- stats
	}
	out
}
