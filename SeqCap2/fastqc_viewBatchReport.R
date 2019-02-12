#' Show fastQC summary results for multiple samples in a 
#' concise report

#' Function is defined first and code is run at the very bottom.

# --------------------------------------------------------------
# Helper routines

#' process basequal
#' generates the plot for base-level quality
processBaseQual <- function(dat,sampName) {
xlims <- c(1,nrow(dat))
plot(0,0,xlim=xlims, ylim=c(0,36),
type='n',xlab="base #",ylab="quality",
xaxt="n",yaxt="n",bty='n',las=1)

# formatting like fastqc
axis(1,at=1:nrow(dat),labels=dat[,1],cex.axis=0.8)
axis(2,at=1:40,labels=1:40,cex.axis=0.8,las=1)
rect(xlims[1],0,xlims[2],20,col='lightpink',border=NA)
rect(xlims[1],20,xlims[2],28,col='lightgoldenrod',border=NA)
rect(xlims[1],28,xlims[2],40,col='darkseagreen1',border=NA)

# plot signal like a boxplot
rect(xleft=(1:nrow(dat))-0.2, xright=(1:nrow(dat))+0.2, 
ybottom=dat[,"Lower Quartile"],ytop=dat[,"Upper Quartile"],
col='yellow')
segments(x0=(1:nrow(dat))-0.2,x1=(1:nrow(dat))+0.2,
y0=dat$Median,lwd=2,col='red')
segments(x0=1:nrow(dat), y0=dat[,"10th Percentile"],
y1=dat[,"90th Percentile"])
title(sampName)
}

# --------------------------------------------------------------
# Work begins


#' parse base quality 
#'
#' @param inDir (char) directory with *_fastqc folders. Expects paired-end
#' reads with "_R1_" and "_R2_" as part of the names
#' @return Side effect of plotting all base quality estimates to a pdf file
fastqc_getBatchReport <- function(inDir) {
	dirList <- dir(inDir,pattern="_fastqc$") 
	cat(sprintf("Got %i files\n",length(dirList)))

		fastQC_sections <- c(
			">>Basic Statistics",
			">>Per base sequence quality",
			">>Per tile sequence quality",
			">>Per sequence quality scores",
			">>Per base sequence content",
			">>Per sequence GC content",
			">>Per base N content",
			">>Sequence Length Distribution",
			">>Sequence Duplication Levels",
			">>Overrepresented sequences",
			">>Adapter Content",	
			">>Kmer Content")
	
	# pass/fail various tests
	testStatus <- matrix("",nrow=length(dirList),ncol=length(fastQC_sections))
	rownames(testStatus) <- dirList
	colnames(testStatus) <- sub(">>","",fastQC_sections)

	pdf(sprintf("%s/fastQC_report.pdf", inDir),width=11,height=11)
	par(mfrow=c(3,2))
	tryCatch({

	ctr <- 1
	for (curd in dirList) {		
		cat(sprintf("%s\n",curd))
		dat <- scan(what="character",
			file=sprintf("%s/%s/fastqc_data.txt",inDir,curd),
			sep="\n",quiet=TRUE)
		idx <- grep(">>Per base sequence quality",dat)
		idx2 <- which(dat==">>END_MODULE")[2]
		basequal <- dat[(idx+1):(idx2-1)]
		# convert to table
		x <- strsplit(basequal,"\t")
		df <- do.call("rbind",x[-1])
		colnames(df) <- x[[1]]
		df <- as.data.frame(df,stringsAsFactors=F)
		for (k in 2:ncol(df)) df[,k] <- as.numeric(df[,k])
		processBaseQual(df,curd)

		for (curt in 1:length(fastQC_sections)) {
			idx <- grep(fastQC_sections[curt],dat)
			x <- strsplit(dat[idx],"\t")[[1]][2]
			testStatus[ctr,curt] <- x
		}
		ctr <- ctr+1
	}
	},error=function(ex) {
		print(ex)
	},finally={
		dev.off()
	})

	# plot test status
	status2 <- matrix(nrow=nrow(testStatus),ncol=ncol(testStatus))
	status2[which(testStatus=="pass")] <- 1
	status2[which(testStatus=="warn")] <- 2
	status2[which(testStatus=="fail")] <- 3
	require(plotrix)
	pdf(sprintf("%s/testStatus.pdf",inDir),height=11,width=6)
	par(mar=c(1,6,7,1))
	tryCatch({
		 plotrix::color2D.matplot(status2,show.values=FALSE,axes=FALSE,
			cs1=c(0,1,1),
			cs2=c(1,1,0),
			cs3=c(0,0,0), xlab="",ylab="")
      axis(3,at=(1:ncol(testStatus))-0.5,
			labels=sub("[Ss]equence","seq",
					sub("Overrepresented","Overrep.",
						sub("[Ll]ength","len.",
						sub("[Qq]uality","qual.",
							colnames(testStatus))))),
            tick=F,cex.axis=0.8,line=-1,las=2)                      
      axis(2,at=seq_len(nrow(testStatus))-0.5,                            
      	labels=sub("_001_fastqc","",rev(rownames(testStatus))),tick=F, 
        las=1,cex.axis=0.7)    
	},error=function(ex) {
		print(ex)
	},finally={
		dev.off()
	})

}

inDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/fastqc"
fastqc_getBatchReport(inDir)
