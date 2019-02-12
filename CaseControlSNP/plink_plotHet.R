
# plot mean heterozygosity per individual
# write ID of those with heterozygosity exceeeding 3SD from mean

args	<- commandArgs(TRUE)
hFile	<- args[1]
pdfFile	<- sprintf("%s.pdf",hFile)

het	<- read.table(hFile,header=T)
het$meanHet <- (het$N.NM. - het$O.HOM.)/het$N.NM. 
pdf(pdfFile)

tryCatch({
	hist(het$meanHet, col="green", 
		 ylab="Heterozygosity Rate", xlab="Heterozygosity", 
		 main= "Distribution of heterozygosity per individual")
	#red lines  +/- 2 standard deviations from the mean
	linea<- mean(het$meanHet) + 2*sd(het$meanHet)
	lineb<- mean(het$meanHet) - 2*sd(het$meanHet)
	abline(v=linea, col="red") 
	abline(v=lineb, col="red")

	#blue lines  +/-3 standard deviations from the mean
	linea<- mean(het$meanHet) + 3*sd(het$meanHet)
	lineb<- mean(het$meanHet) - 3*sd(het$meanHet)
	abline(v=linea, col="blue") 
	abline(v=lineb, col="blue")
}, finally={
	dev.off()
})

# write list of those failing to file.
errora<- subset(het, het$meanHet > mean(het$meanHet) + 3*sd(het$meanHet))
errorb<- subset(het, het$meanHet < mean(het$meanHet) - 3*sd(het$meanHet)) 
errorall<-rbind(errora[1:2], errorb[1:2])
outFile <- sprintf("%s/fail-het-qc.txt", dirname(hFile))
write.table(errorall, file=outFile, col.names=FALSE, 
			row.names=FALSE, quote=FALSE)




