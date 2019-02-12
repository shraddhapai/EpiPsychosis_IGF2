
###inFile <- "~/Epigenetics/NARSAD/anno/MGI_Geno_DiseaseDO_C57BL6.txt"
###dat <- read.delim(inFile,sep="\t",h=F,as.is=T)
###colnames(dat)[8] <- "DO.Disease.ID"
###colnames(dat)[9] <- "OMIM.ID"

inFile <- "~/Epigenetics/NARSAD/anno/mouse/MGI_DO.rpt"
mgiFile <- "~/Epigenetics/NARSAD/anno/mouse/MRK_List2.rpt"
dat <- read.delim(inFile,sep="\t",h=T,as.is=T)
dat <- subset(dat,Mouse.MGI.ID!="")

mgi <- read.delim(mgiFile,sep="\t",h=T,as.is=T)

outFile <- sub("rpt","gmt",basename(inFile))
system(sprintf("cat /dev/null >%s", outFile))
pathList <- list()
ctr <- 0
for (k in unique(dat$OMIM.ID)) {
	idx <- which(dat$OMIM.ID %in% k)
	gn <- unique(dat$Mouse.MGI.ID[idx])
	gn2 <- mgi$Marker.Symbol[which(mgi$MGI.Accession.ID %in% gn)]
	cat(sprintf("%s: %i entries\n",k,length(gn2)))
	if (length(gn2)>=10 & length(gn2) <=200) {
		ctr <- ctr+1
		cat(sprintf("%s\t%s\t%s\n",k,k,paste(gn2,collapse="\t")),
			file=outFile,append=TRUE)
	}
}
cat(sprintf("%i pathways between 10 and 200\n",ctr))



