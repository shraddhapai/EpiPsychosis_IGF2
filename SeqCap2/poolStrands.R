#' pools reads from + and - strands and assigns * as strand. 
#' @param recFile (char) path to *records.Rdata which contains base-level
#' methylation on a per-target level.
#' @param getTargetGR (logical) if TRUE, returns target_gr as output
#' @return if getTargetGR=FALSE (list) one key per target GRanges. 
#'		chr,pos,strand,C_count,CT_count.  Counts are pooled across strands.
#' if getTargetGR=TRUE (list of length 2) 1=strand-pooled methylation list as
#' above. 2=target_GR
poolStrands <- function(recFile,getTargetGR=FALSE)  {
	cat("Pooling strands\n")
	lnames <- load(recFile)
	out <- list()
	for (k in 1:length(rec)) {
	cur <- rec[[k]]
	if (!is.null(dim(cur))) {
		cur <- cur[,c(1:3,7:8)]
		cur$strand <- as.character(cur$strand)
		cur$pos <- as.integer(as.character(cur$pos))
		midx <- which(cur$strand=="-")
		cur$strand <- "*"
		if (any(midx)) {
			cur$pos[midx] <- cur$pos[midx]-1 # shift to + coord
			cur <- aggregate(cur[,c("C_count","CT_count")],
				by=list(chr=cur$chr,pos=cur$pos,strand=cur$strand),
				FUN=sum)
		}
	}
	out[[k]] <- cur
	}
	names(out) <- names(rec)
	if (getTargetGR) { return(list(rec=out,target_GR=target_gr)) 
	} else { return(out)
	}
}

