# plot missing individual data output from plink.

args	<- commandArgs(TRUE)
mFile	<- args[1];

pdfFile	<- sprintf("%s.pdf", mFile)
IMISS	<- read.table(mFile,header=T,as.is=T)
pdf(pdfFile);
tryCatch({
	plot( (1:dim(IMISS)[1])/(dim(IMISS)[1]-1), sort(1-IMISS$F_MISS), 
		 main="Individual call rate cum. distribution",
		 xlab="Quantile", ylab="Call Rate" ); 
	grid()
	abline(h=0.999,lty=3,lwd=2)
},error=function(ex) { 
	print(ex)
}, finally={
	dev.off()
})
