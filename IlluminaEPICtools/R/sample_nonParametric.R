#' get a pvalue by sampling with replacement from background
#' 
#' @param fg (numeric) foreground set
#' @param bg (numeric vector) background set
#' @param FUN (function) function to compare
#' @param type (char) [twotailed|fg_greater|fg_less] type of tailed test
#' @param showPlot (logical) show real and null distributions
#' @param numCores (integer) for parallel processing
#' @import foreach
#' @export
sample_nonParametric <- function(fg,bg,R=1000,FUN=median,type="twotailed",
	showPlot=TRUE,numCores=1L) {
	x <- numeric() 
	val <- FUN(fg)
	
	require(doParallel)
	require(foreach)
	cl <- makeCluster(numCores)
	registerDoParallel(cl)

	x  <- foreach (r=1:R) %dopar% {
		FUN(sample(bg,length(fg),replace=TRUE))
	}
	x <- unlist(x)

	stopCluster(cl)

		pval <-	switch(type,
			twotailed=min(sum(x<=val),sum(x>=val)),
			fg_greater=sum(x >= val),
			fg_less=sum(x <= val),
			stop("invalid type"))
	

		a <- min(x,val); b <- max(x,val)
		plot(density(x),xlim=c(a*0.99,b*1.01),
			 main=sprintf("%i Fore; %i Back; R=%i (p=%1.2e)",
						  length(fg),length(bg),R,pval/R)); 
		abline(v=val,lwd=3,col='red')
		if (pval< 1) pval <- 1/R else pval  <- pval/R
		pval
}

