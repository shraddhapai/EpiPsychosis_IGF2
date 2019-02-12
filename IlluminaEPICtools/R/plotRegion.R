#' plot methylation for a region
#' 
#' @param mset (MethylSet-Genomic)
#' @param rgn (GenomicRanges) range to show
#' @param gp (factor) value for each sample
#' @param buffer (integer) +/-bp to plot around rgn
#' @param showGenes (logical) show gene bodies?
#' @param gene_GR (Granges) ranges of genes, if showGenes=TRUE
#' @export
#' @import RColorBrewer
plotRegion <- function(mset,rgn,gp,buffer=1000,showGenes=TRUE,
					   gene_GR) {

	orig_rgn <- rgn
	rgn <- resize(rgn, width=width(rgn)+(2*buffer),fix="center")

	x1 <- start(rgn)
	x2 <- end(rgn)

	gp <- as.factor(gp)
	locs <- getLocations(mset)
	ol <- findOverlaps(locs,rgn)
	idx <- unique(ol@queryHits)
	betas <- getBeta(mset[idx,])

	pal <- brewer.pal(n=length(levels(gp))+1, name="Spectral")
	pal <- pal[-1]

	
	ylims <- c(min(betas),max(betas))
	plot(NA,NA,xlim=c(x1,x2),
		 ylim=c(-0.1,1),las=1,bty='n',
		 xlab="Position(bp)", ylab="Beta",
		 mar=c(4,2,4,2))
	xpos <- start(locs)[idx]
	for (k in 1:ncol(betas)) {
		points(xpos,betas[,k],col=pal[gp[k]],type="p",pch=16,cex=0.3)		
	}

	# add a group trendline
	ctr <- 1
	for (g in levels(gp)) {
		mu <- rowMeans(betas[,which(gp%in%g),drop=FALSE])
		ks <- ksmooth(xpos,mu,"normal",bandwidth=max(10,0.01*width(rgn)))
		lines(ks$x,ks$y,col=pal[ctr])
		ctr <- ctr+1
	}
	
	if (showGenes) {
		gn <- subsetByOverlaps(gene_GR,rgn)
		rect(xleft=start(gn),xright=end(gn),ybottom=-0.1,ytop=0,
			 col='orange')
		starts <- start(gn)
		str <- as.character(strand(gn))
		pidx <- which(str%in%"+")
		midx <- which(str%in%"-")
		alen <- max(10,0.01*width(rgn))
		if (any(pidx)) {
			arrows(x0=starts[pidx],x1=starts[pidx]+alen,
				y0=0.05,y1=0.05,code=2)
		}
		if (any(midx)) {
			starts[midx] <- end(gn)[midx]
			arrows(x0=starts[midx],x1=starts[midx]-alen,
				   y0=0.05,y1=0.05,code=2)
		}
		segments(x0=starts,x1=starts,y0=0,y1=0.05)

	}
	rect(xleft=start(orig_rgn),xright=end(orig_rgn),
		 ybot=-0.1,ytop=0,col='green')

	legend("topright",legend=levels(gp),fill=NA,pch=16,col=pal,bty='n',
		   border=NA)
	title(sprintf("%s: %i-%i (%i kb)", seqnames(rgn),start(rgn),end(rgn),
				  round(width(rgn)/1000)))
	cat(sprintf("%s:%i-%i",seqnames(rgn),start(rgn),end(rgn)))

}
