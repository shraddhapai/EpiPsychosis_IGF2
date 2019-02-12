# plot MDS of data genotypes against HapMapv3 genotypes
rm(list=ls())
require(combinat)

loc <- "LOCAL" # LOCAL | SCINET

# ----------------------------------
# Dataset specific files
args	<- commandArgs(TRUE)
dataFam <- args[1]
mdsFile <- args[2]
###HMdir <- args[3]

###if (loc == "SCINET") {
###	# .fam file for our data after cleaning
###	dataFam	<- "/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399/hg19/SUS19399.hg19.sorted" # args[1];
###	# output of plink --mds call
###	mdsFile	<- "/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399/hg19/MDS.mds" #args[2];
###	# dir with HapMap3 fam files
	HMdir	<- "/home/shraddhapai/imputation_annotation/hapmap3_pop"
####mdsFile	<- "/home/spai/BaderLab/PNC/plinkQC/GO_v3/MDS.mds"
###} else {
###	dataFam <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlSNP/plinkQC/SUS19399.hg19.sorted.fam"
####	mdsFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlSNP/plinkQC/MDS.mds"
###	mdsFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlSNP/plinkQC/HMDATA_PCA.eigenvec"
###	HMdir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/hapmap3_pop"
###}



# blue colours for recent african origin
# warm colours for east asian origin
# brown for indian 
# purple/pink for european origin
# green for mesoamerican
clrMap	<- list( 
  ASW="#74a9cf", # african amer in texas
  YRI="#023858", # yoruban 
  LWK="#02818a", # luhya, kenyan
  MKK="#016c59", # Masai, Kenya
  CHB="#bd0026", # han chinese
  CHD="#ffeda0", # chinese in metro LA
  JPT="#feb24c", # japanese in tokyo
  GIH="#993404", # gujarati indian in texas
  MEX="#fdae6b", # mexican
  TSI="#dd3497", # toscani italian
  CEU="#762a83" # European descent, utah
 )
HM_POP	<- c("ASW","CHB","GIH","JPT","MEX","TSI","CEU","CHD","LWK","MKK","YRI")


# ----------------------------------
# Work starts here.
outDir <- sprintf("%s/ethnicities", dirname(dataFam))
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

outFile	<- sprintf("%s/fail-divergent-ancestry.txt", outDir)
ethFile	<- sprintf("%s/HapMap3_based_ancestry.txt",outDir)
system(sprintf("cat /dev/null > %s", outFile))

if (any(grep(".mds$",mdsFile))) {
	fileType <- "MDS" 
	mds	<- read.table(mdsFile,header=TRUE)
} else {
	fileType <- "PCA"
	mds <- read.table(mdsFile,h=F,as.is=T)
	colnames(mds)[1:ncol(mds)] <- c("FID","IID",
									paste("C",1:(ncol(mds)-2),sep=""))
}
cat(sprintf("Reduction type provided: %s\n", fileType))

dat	<- read.table(dataFam)

samp <- merge(dat,mds, by.x=c("V1","V2"),by.y=c("FID","IID"))
samp$samp_eth <- NA
n	<- length(HM_POP)

cat("* Read in population .fam files\n")
popList <- list()
pop	<- character()

outAll <- sprintf("%s/MDS.ethnicity.all.txt",outDir)
system(sprintf("cat /dev/null > %s", outAll))

sclr <- "red"

pdfFile <- sprintf("%s.pdf", mdsFile)
pdf(pdfFile)
#par(mfrow=c(2,2))
tryCatch({
comps <- combn(5,2) # plot first five PCs
	# plot reference populations
	for (combNum in 1:ncol(comps)) {
		cx <- paste("C",comps[1,combNum],sep="")
		cy <- paste("C",comps[2,combNum],sep="")

		for (p in 1:length(HM_POP)) {
			if (combNum==1) cat(sprintf("%s\n",HM_POP[p]))
			x	<- read.table(sprintf("%s/%s.fam", HMdir, HM_POP[p]))
			y	<- merge(x, mds, by.x=c("V1","V2"),by.y=c("FID","IID"))
			clr	<- clrMap[[HM_POP[p]]];
			if (p==1) {
				plot(y[,cx],y[,cy], xlab=cx,ylab=cy, 
					 xlim=c(min(mds[,cx]),max(mds[,cx])),
					 ylim=c(min(mds[,cy]),max(mds[,cy])),
					 col=clr,cex=0.6,bty='n',
					 main=sprintf("%s\n %s vs %s from HapMap3 genotypes",
								  fileType,cy,cx,basename(fileType)))
	#	rect(lim[1],lim[3],lim[2],lim[4],border=clrMap[[nm]],lty=3,
#			 lwd=2)
#		limclr <- c(limclr, clrMap[[nm]])
			} else {
				points(y[,cx],y[,cy],col=clr,cex=0.6)
			}
			
			#points(samp[,cx],samp[,cy],col=sclr,cex=0.7,pch=4)
			text(samp[,cx],samp[,cy],labels=samp$V2,col='grey40',
				 cex=0.5)
			pop <- c(pop,clr)
			popList[[HM_POP[p]]] <- y
	
			# get threshold for pc1 and pc2
			mysd <- c(sd(y[,8]), sd(y[,9]));
			mymu <- c(mean(y[,8]),mean(y[,9]));
			tmp	<- c( mymu-(3*mysd), mymu+(3*mysd))
			lim <- tmp[c(1,3,2,4)]
			# assign sample ancestry if within bounds of this population
			if (combNum==1) {
				pc1	<- samp[,cx]>lim[1] & samp[,cx] < lim[2]
				pc2 <- samp[,cy]>lim[3] & samp[,cy] < lim[4]
				idx	<- which(pc1 & pc2)
				if (any(idx)) {
					samp$samp_eth[idx] <- HM_POP[p]
					cat("********\n")
					print(table(samp$samp_eth,useNA="always"))
					cat("********\n")
				#text(samp[idx,cx],samp[idx,cy],samp[idx,2],
				#cex=0.7,col='blue')
			
				cat(sprintf("\t%s: %i individuals\n", HM_POP[p], length(idx)))
				write.table(samp[idx,1:2],
					file=sprintf("%s/%s.txt", outDir,HM_POP[p]), 
					col=F,row=F,quote=F)
			
				write.table(cbind(samp[idx,1:2],HM_POP[p]),file=outAll,
					sep="\t",col=F,row=F,quote=F,append=TRUE)
				}
			legend("topleft", pch=c(rep(20,n),4),
				col=c(pop, 'red'),ncol=2,legend=c(HM_POP,"study"),
				cex=0.7)
			print(table(samp$samp_eth,useNA="always"))
			}
		}

		if (combNum==1) {
		# mark samples not in any group, divergent
		na_idx	<- which(is.na(samp$samp_eth))
		cat(sprintf("\tOutliers: %i individuals\n", length(na_idx)))
		if (any(na_idx)) {
			#points(samp[na_idx,cx],samp[na_idx,cy],col=sclr,
			#   cex=0.4,pch=4)
			write.table(samp[na_idx,1:2], file=outFile,col=F,row=F,quote=F)
        	write.table(cbind(samp[na_idx,1:2],"divergent"),file=outAll,
               sep="\t",col=F,row=F,quote=F,append=TRUE)

		tmp	<- samp[,c("V1","V2","samp_eth")]
		write.table(tmp,file=ethFile,col=FALSE,row=FALSE,quote=FALSE)
		print(table(samp$samp_eth,useNA="always"))
		}}
	}
		

	# legend

},error=function(ex){
	print(ex)
},finally={
	dev.off()
})

cat(sprintf("Outliers written to:\n%s\n", outFile))
cat(sprintf("Ethnicity written to:\n%s\n", ethFile))
