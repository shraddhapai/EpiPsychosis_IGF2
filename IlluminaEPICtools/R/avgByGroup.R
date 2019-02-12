#' average beta values by group
#'
#' @param b (matrix) beta values
#' @param g (factor) grouping variable
#' @return (matrix) mean of values belonging to same group
#' @export
avgByGroup <- function(b,g) {
	lv <- levels(g)
	b2 <- matrix(NA,nrow=nrow(b),ncol=length(lv))
	ctr <- 1
	for (nm in lv) {
		b2[,ctr] <- rowMeans(b[,which(g %in% nm),drop=FALSE])
		ctr <- ctr+1
	}
	rownames(b2) <- rownames(b)
	colnames(b2) <- lv

	b2
}
