# plot IBS and generate list of indiv with IBS > threshold.

args <- commandArgs(TRUE)
ibsFile	<- args[1]
#ibsFile <- "/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399_New/hg19/SUS19399.hg19.sorted_IBS.genome"
print(ibsFile)

cat("* Reading IBS file. This could take a while\n")
IBS		<-read.table(ibsFile, h=T)
pdfFile	<- sprintf("%s.pdf", ibsFile)
pdf(pdfFile)
tryCatch({
	hist(IBS$PI_HAT, col="green", breaks=100,  
		 xlab= "Estimated mean pairwise IBD", 
		 main="Pairwise Identity-by-descent"
	)
	#related individuals are assumed to have a pairwise IBD of over 0.185
	abline(v=c(0.185)) 
}, error=function(ex){
	print(ex)
}, finally={
	dev.off()
})
# Still in R. Make list
out<- subset(IBS, IBS$PI_HAT>0.185) #pick first individual in pair
write.table(out[,c(1,2)], file=sprintf("%s/fail-IBD.txt",dirname(ibsFile)),
			col.names=F, row.names=F, quote=F)

