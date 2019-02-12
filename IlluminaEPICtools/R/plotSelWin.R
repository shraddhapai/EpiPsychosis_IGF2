#' plots a signal in windows of interest
#'
#' @param gr (GRanges) signal of interest
#' @param cname (char) column with signal to plot
#' @param targetWin (GRanges) windows in which to plot the signal
#' @param refLine (numeric) yintercept for horizontal reference
#' line
#' @param ... params for plot()
#' @export
plotSelWin <- function(gr, cname,targetWin,refLine=0,...) {
	if ( length(targetWin)> 1)par(mfrow=c(3,1),las=3)
	for (i in 1:length(targetWin)){
		win <- targetWin[i]
		cur <- subsetByOverlaps(gr,win)
		df <- elementMetadata(cur)
		plot(start(cur)/1000,df[,cname],type="o",
			 xlim=c(start(win)/1000,end(win)/1000),
			 xlab="Coordinate (kb)",ylab=cname,bty='n',...)
		points(start(cur)/1000,df[,cname],col='red',pch=16,
			   cex=0.4)
		title(sprintf("%s:%s-%s (%1.2f kb)",
			as.character(seqnames(win)),
			prettyNum(start(win),big.mark=","),
			prettyNum(end(win),big.mark=","),
			width(win)/1000),cex=1.5)
		abline(h=refLine,lty=3,col='green',lwd=2)
	}
}
