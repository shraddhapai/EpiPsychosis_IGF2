# plot locus-level methylation of neun+ vs neun-

indir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/methylation/cph/locuswise"

fList <- dir(indir,pattern="Rdata")
minCvg <- 10
out <- list()
for (fName in fList) {
	print(fName)
	load(sprintf("%s/%s",indir,fName))
	fTemp <- sub("records.Rdata","_cph.bed",fName)
	system(sprintf("cat /dev/null > %s", fTemp))
	meth <- lapply(rec, function(x) {
		if (is.na(x)) return(NULL)
		x <- subset(x,CT_count >= minCvg)
		if (nrow(x) > 0) {
			pctm <- round((x$C_count/x$CT_count)*100,digits=2)
			tmp <- cbind(x[,c(1,2,2)],".",pctm,x[,3])
			write.table(tmp,file=fTemp,sep="\t",col=F,row=F,quote=F,
				append=TRUE)
			#return(tmp)
			return(sum(x$C_count)/sum(x$CT_count))
		} else {
			return(NULL)
		}
	})
		#write.table(out[[fName]],file=fTemp,sep="\t",col=F,row=F,quote=F)
		#out[[fName]] <- meth
}
