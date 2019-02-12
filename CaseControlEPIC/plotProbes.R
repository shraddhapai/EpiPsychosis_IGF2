#' plot sequence of probes and highlight samples with max min
#' @param betas
#' @param pd 
#' @param outFile
plotProbes_showSamp <- function(betas, pd,outFile) {
	x <- betas[,which(pd$DX %in% "case")]
	y <- betas[,which(pd$DX %in% "control")]
	x2 <- sort(colMeans(x))
	y2 <- sort(colMeans(y))

	nr <- nrow(betas)
		
	pdfFile <- sub(".txt",".pdf",outFile)
	pdf(pdfFile,width=11,height=6)
	tryCatch({
	plot(0,0,type='n',xlim=c(1-0.2,nr+0.2),ylim=c(0,1),
		 ylab="Beta values",las=1,bty='n')
	title(sprintf("Sample-level beta view\n{ %s }", 
		  paste(rownames(betas),collapse=",")))
	legend("topleft", legend=c("case, all", "control, all",
		"case,bottom","control, top"), fill=NA,border=NA,
		   col=c("pink","lightskyblue2","red","blue"),
		   pch=c(1,4,1,4),cex=0.8,lwd=1,bty='n')

	# plot all cases
	for (k in 1:ncol(x)) {
		points(1:nr,x[,k],type='o',col=rgb(1,0.7,0.7,0.5))
	}
	# plot all controls
	for (k in 1:ncol(y)) {
		points(1:nr,y[,k],type='o',col=rgb(0.7,0.7,1,0.5),pch=4,cex=0.8)
	}
	# lowest of the cases
	for (k in 1:6) {
		points(1:nr,x[,names(x2)[k]],type='o',col='red')
	}
	cat("-----------\n")
	idx <- which(rownames(pd)%in% names(x2)[1:6])
	cat(sprintf("Lowest cases: { %s }\n", 
				paste(pd$Sample_Name[idx],collapse=",")))
	print(table(pd$DX[idx]))
	# highest of the controls
	ns <- length(y2)
	for (k in 0:5) {
		points(1:nr,y[,names(y2)[ns-k]],type='o',
			   pch=4,cex=0.8,col='darkblue')
	}
	}, error=function(ex){
		print(ex)
	}, finally={
		dev.off()
	})

	cat("-----------\n")
	idx <- which(rownames(pd)%in% names(y2)[ns+(-5:0)])
	cat(sprintf("Highest controls: { %s }\n", 
				paste(pd$Sample_Name[idx],collapse=",")))
	print(table(pd$DX[idx]))

	cat("-----------\n")
	idx <- which(rownames(pd)%in% names(x2)[ns+(-5:0)])
	cat(sprintf("\nHighest cases: { %s }\n", 
				paste(pd$Sample_Name[idx],collapse=",")))
	print(table(pd$DX[idx]))

	cat("-----------\n")
	idx <- which(rownames(pd)%in% names(y2)[1:6])
	cat(sprintf("Lowest controls: { %s }\n", 
				paste(pd$Sample_Name[idx],collapse=",")))
	print(table(pd$DX[idx]))

	system(sprintf("cat /dev/null > %s",outFile))
	cat("Cases\n",file=outFile,append=TRUE)
	tmp <- pd[which(pd$DX %in% "case"),]
	midx <- match(rownames(tmp), colnames(x))
	x <- x[,midx]
	colnames(x) <- tmp$Sample_Name
	suppressWarnings(write.table(t(x),file=outFile,sep='\t',
			col=T,row=T,quote=F,append=TRUE))

	cat("Controls\n",file=outFile,append=TRUE)
	tmp <- pd[which(pd$DX %in% "control"),]
	midx <- match(rownames(tmp), colnames(y))
	y <- y[,midx]
	colnames(y) <- tmp$Sample_Name
	suppressWarnings(write.table(t(y),file=outFile,sep='\t',
			col=T,row=T,quote=F,append=TRUE))

}
