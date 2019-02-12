#' plot betas in window around selected probes
#' 
#' @param mset (MethylSet)
#' @param selProbes (char) probes to plot
#' @param winWidth (numeric) window width centered on selProbes
#' @param writeMelted (logical) if TRUE write melted matrix to file.
#' has probe name and sample info
#' @param outPref (char) prefix for output file
#' @export
#' @import reshape2
dmp_plotLocalBeta <- function(mset,selProbes,group=NULL,winWidth=20000,
							  writeMelted=FALSE,outPref=NULL){
	locs <- getLocations(mset)
	pd <- pData(mset)
	pd$NAME <- rownames(pd)
	pd <- pd[,c("NAME","SEX","DX","DIST.DX")]
	for (cur in selProbes) {
		print(cur)
		idx <- which(names(locs) %in% cur)
		
		target_win <- resize(locs[cur],fix="center",width=winWidth)
		win <- subsetByOverlaps(locs,target_win)
		cat("getting distance\n")
		d <- as.data.frame(distanceToNearest(win,locs[cur]))
		cat(sprintf("got %i rows\n",nrow(d)))
		d <- subset(d, distance <= floor(winWidth/2))
		cat(sprintf("filtered = %i rows\n", nrow(d)))

		b <- getBeta(mset[names(win)])
		m <- melt(b)
		colnames(m) <- c("PROBE_NAME","NAME","VALUE")
		m[,1] <- as.character(m[,1])

		# merge betas with pheno table
		x <- merge(pd, m, by="NAME")
		cur_locs <- as.data.frame(locs[names(win)])
		cur_locs$PROBE_NAME <- rownames(cur_locs)
		y <- merge(x=x,y=cur_locs,by="PROBE_NAME")

		if (writeMelted) {
			if (is.null(outPref)) stop("outPref is NULL")
			outFile <- sprintf("%s.meltedBeta_%s.txt",outPref,cur)
			options(scipen=10)
			write.table(y,file=outFile,sep="\t",col=TRUE,row=FALSE,quote=F)
		}

		# use ggplot to plot beta over region, with a smoothed line
		# condition on case/control, male/female, ancestry etc.,
		require(ggplot2)
		y$start <- y$start/1000
		y <- as.data.frame(y)
		p <- ggplot(y, aes(start, VALUE))
		p <- p + stat_smooth(aes(colour=DX)) + 
			geom_point(aes(colour=DX),cex=0.3)
		p <- p + ylab("beta") + xlab("Location (kb)")
		p <- p + ggtitle(sprintf("%s\n%s: %i-%i",cur,
								 as.character(seqnames(target_win)[1]),
								 start(target_win),end(target_win)))
		print(p)

	}
}
