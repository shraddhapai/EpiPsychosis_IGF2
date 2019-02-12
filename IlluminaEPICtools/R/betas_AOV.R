#' run ANOVA across groups
#'
#' @param betas (matrix) probe-wise beta values
#' @param groupv (factor) group values, corresponding to columns of betas
#' @param topX (integer) limit to top variable probes. If NULL, uses all
#' @param numCores (integer) for parallel processing
#' @return (list) topVar: probe IDs for top variable probes
#' @import doParallel
#' @import foreach
#' @export
betas_AOV <- function(betas, groupv,topX=10000L,numCores=8L) {
	require(doParallel)
	cl <- makeCluster(numCores)
	registerDoParallel(cl)
	print(dim(betas)) # namespace leak?

	if (!is.null(topX)) {
		cat("Computing probe-wise variability\n")
		print(system.time(bsd <- foreach(k=1:nrow(betas)) %dopar% {
			sd(betas[k,])	
		}))
		cat(sprintf("\tLimiting to top %i\n",topX))
		idx <- order(bsd,decreasing=TRUE)[1:topX]
		betas <- betas[idx,]
	}

	browser()

	out <- foreach(k=1:nrow(betas)) %do% {
		fit <- aov(betas[k,]~groupv)
		browser()
	}

}

