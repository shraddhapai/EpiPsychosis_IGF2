#' run FET to ascertain overep. of sig probes
#'
#' @param fg_GR (GRanges) foreground
#' @param bg_GR (GRanges) background
#' @param targets (data.frame or list) if char, path to bed files 
#' containing genomic regions (Type, description, File); 
#' if list assumes a list of GRanges
#' @import GenomicRanges
#' @return (matrix) one row per data type. num_fg, num_bg, pct_fg, pct_bg, pvalue from FET
#' @export
dmp_Fisher <- function(fg_GR, bg_GR, targets) {
	if (class(targets)=="data.frame") {
		outmat <- matrix(nrow=nrow(targets),ncol=5)
		rownames(outmat) <- targets[,1]
	} else {
		outmat <- matrix(nrow=length(targets),ncol=5)
		rownames(outmat) <- names(targets)
	}
	colnames(outmat) <- c("num_fg","num_bg","pct_fg","pct_bg","FET_p")
	ctr <- 1
	
	if (class(targets)=="data.frame"){
		loopOver <- targets$File
	} else {
		loopOver <- names(targets)
	}

	for (f in loopOver) {
		if (class(targets) %in% "data.frame") {
			cat(sprintf("\t%s\n",targets[ctr,1]))
			bed <- read.delim(f,sep="\t",h=F,as.is=T)
			bed_GR <- GRanges(bed[,1],IRanges(bed[,2],bed[,3]))
		} else {
			cat(sprintf("\t%s\n", f))
			bed_GR <- targets[[f]]
		}
		seqlevels(bed_GR, force=TRUE) <- seqlevels(fg_GR)

		fg_ol <- subsetByOverlaps(fg_GR, bed_GR)
		bg_ol <- subsetByOverlaps(bg_GR, bed_GR)
		
		nF <- length(fg_ol); nB <- length(bg_ol)
		fmat <- matrix(nrow=2,ncol=2)
		fmat[1,1]<-nF  					# FG, yes
		fmat[1,2]<-length(fg_GR)-nF    	# FG, no
		fmat[2,1]<-nB    				# BG, yes
		fmat[2,2]<-length(bg_GR)-nB    	# BG, no
		fet <- fisher.test(fmat)

		outmat[ctr,] <- c(nF,nB,nF/length(fg_GR),nB/length(bg_GR),
						  fet$p.value)
		ctr <- ctr+1
	}
	
	# plot % fg, % bg, and asterisks showing overrep
	cat("* Making barplot\n")
	dat <- outmat[,3:4]*100
	par(las=1,cex.axis=1)
	x <- barplot(t(dat),beside=T,col=c("red","grey50"),
				 border='white',ylim=c(0,max(dat)*1.3),
				 ylab="Proportion (%)")
	title("% foreground/background probes: FET results")
	x <- colMeans(x)

	text(x,max(dat)*1.2,sprintf("%1.2e",outmat[,5]),font=2)
	idx <- which(outmat[,5] < 0.001)
	if (any(idx)) text(x[idx],max(dat)*1.1,"***",col='red',font=2)
		
	outmat
}
