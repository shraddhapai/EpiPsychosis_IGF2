#' compute row-wise variance
#'
#' @param b (matrix)
#' @export
getVariance <- function(b) {
	mu	<- rowMeans(b)
	n	<- ncol(b)
	# corrected sd estimator
	x	<- rowSums((b-mu)^2)/(n-1)
	x
}
