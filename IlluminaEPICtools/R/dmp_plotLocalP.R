#' plot pvalue around hits, generate significant region list
#'
#' @details A window of minWin is examined around each significant 
#' probe (selProbes). If it contains probes with FDR < FDR2, the
#' window is extended to the maximum distance of such a probe,
#' with a buffer of 10%.
#' If no such probes are present, the window has a minimum width
#' of minWin.
#' @param mset (MethylSet)
#' @param dmpRes (data.frame) named data.frame with 'pval' column
#' @param selProbes (char) probes to plot
#' @param winWidth (numeric) window width centered on selProbes
#' @param FDR1 (num 0-1) probes with higher significance, marked
#' in red in the plot
#' @param FDR2 (num 0-1) larger than FDR1, probes with moderate
#' significance, marked in orange in the plot. Also used to 
#' win
#' @param minWin (integer) the minimum window width around a 
#' signal probe if no other moderate signal probes are present
#' @param buffer (integer) width by which to extend a region
#' @param plotMe (logical) if FALSE does not plot windows
#' @export
#' @return (data.frame) windows around significant probes
dmp_plotLocalP <- function(mset, dmpRes,selProbes,winWidth=20000,
						   FDR1=0.1,FDR2=0.3,minWin=5000,buffer=50,
						   verbose=FALSE){
	locs <- getLocations(mset)
	locs <- locs[which(names(locs)%in% dmpRes$NAME)]
	cat(sprintf("filtered to %i locs\n",length(locs)))

	par(mfrow=c(3,1),las=1,font.axis=2)

	# chrom,start,end, 
	out <- data.frame(seqnames=character(0),
					  start=integer(0),end=integer(0),
					  FDR1=character(0),FDR2=character(0))

	for (cur in selProbes) {
		print(cur)
		idx <- which(names(locs) %in% cur)
		
		target_win <- resize(locs[cur],fix="center",width=winWidth)
		win <- subsetByOverlaps(locs,target_win)
		if (verbose) cat("getting distance\n")
		d <- as.data.frame(distanceToNearest(win,locs[cur]))
		if (verbose) cat(sprintf("got %i rows\n",nrow(d)))
		d <- subset(d, distance <= floor(winWidth/2))
		if (verbose) cat(sprintf("filtered = %i rows\n", nrow(d)))
		wingr <- win[d[,1]]
		win <- as.data.frame(wingr)
		win$NAME <- rownames(win)
		win$distance <- d[,3]
		x <- merge(x=win,y=dmp,by="NAME")
		win <- x[order(x$start),];rm(x)

		sig_probe <- locs[cur]
		ww <-(max(win$end)-min(win$start))/1000
		myq <- win$qval[which(win$NAME %in% cur)]
		txt <- sprintf("%s (Q = %1.2f)\n%s (%i-%i) (%1.2f kb win)",
				cur,myq,as.character(seqnames(sig_probe)),
				start(sig_probe),end(sig_probe),ww)
		
		plot(win$start/1000, -log10(win$pval),
			 type='o',main=txt,lwd=2,pch=16,
			xlab="Probe start (Kb)",ylab="-log10(p)",
			xlim=c(start(target_win)/1000,end(target_win)/1000),
			ylim=c(0,10),cex=0.2)

		idx <- which(win$qval < FDR2)
		win$start <- win$start/1000
		points(win$start[idx],-log10(win$pval)[idx],
			   col='orange',cex=1.3,pch=16)

		idx <- which(win$qval < FDR1)
		points(win$start[idx],-log10(win$pval)[idx],
			   col='red',cex=1.5,pch=16)

		# pick and draw 'signal region'
		sig <- subset(win,qval < FDR2)
		if (nrow(sig)==1) # no other probe in window
			curwin <- resize(locs[cur],fix="center",width=minWin)
		else {
			curwin <- GRanges(sig$seqnames[1],
							  IRanges(min(sig$start*1000),max(sig$end)))
			newwd <- max(minWin, width(curwin)+(buffer*2))
			curwin <- resize(curwin,fix="center",width=newwd)
		}
		rect(xleft=start(curwin)/1000,xright=end(curwin)/1000,
			 ybottom=8,ytop=8.5,col='green',border=NA)

		legend("topleft",legend=c("dmp",
				sprintf("Q<%1.1f",FDR1),sprintf("Q<%1.1f",FDR2)),
			   text.col=c("blue4","red","orange"),
			   col=c("blue4","red","orange"),
			   fill=NA,border=NA,lwd=0,pch=c(18,16,16),
			   bty='n',lty=0,cex=1.3)
		abline(h=-log10(c(0.05,1e-3)),lty=3,col='grey50',lwd=1.5)

		idx <- which(win$NAME %in% cur)
		points(win$start[idx],-log10(win$pval)[idx],pch=18,
			   col="blue4",cex=1.4)

		# add region of interest
		newrow <- data.frame(
					seqnames=as.character(seqnames(curwin)), 
					start=start(curwin),end=end(curwin),
					FDR1=paste(win$NAME[which(win$qval<FDR1)],
							   collapse=","),
						FDR2=paste(win$NAME[which(win$qval<FDR2)],
							   collapse=","))
		for (k in c(1,4,5)) newrow[,k] <- as.character(newrow[,k])
		out <- rbind(out,newrow)
		
	}
	out
}
