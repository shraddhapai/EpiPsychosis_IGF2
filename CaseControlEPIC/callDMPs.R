# ------------------------------------------------------
# call DMPs
# ------------------------------------------------------
rm(list=ls())
require(minfi)
require(IlluminaEPICtools)
source("dmpFinder_DX.R"); 
#source("plotProbes.R"); 

rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398"
cleanDat <- sprintf("%s/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",
	rootDir)
snpPCFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlSNP/plinkQC/SUS19399.hg19.sorted_-clean_PCA.eigenvec"
annoFiles 	<- read.delim("Annotations.txt",sep="\t",h=T,as.is=T)
updatedSampKey <- "/home/shraddhapai/Epigenetics/NARSAD/input_files/NARSAD_sampleKey_171129.txt"

dt <- format(Sys.Date(),"%y%m%d")
outRoot <- rootDir
# ------------------------------------------------------
cat("Loading cleaned data\n")
load(cleanDat)

cat("\n* Preparing for DMP analysis\n")
pd <- pData(MSet.genome)

# update PMI for last three samples
p2 <- read.delim(updatedSampKey,sep="\t",h=T,as.is=T)
pd$PMI <- as.numeric(pd$PMI)
idx <- which(is.na(pd$PMI)) # find missing pmi
cat("Filling in missing PMIs\n")
for (k in idx) {
	idx2 <- which(as.character(p2$Internal.ID) == pd$Sample_ID[k])
	cat(sprintf("\tSample %s\n", pd$Sample_ID[k]))
	pd$PMI[k] <- p2$PMI[idx2]
}
pData(MSet.genome) <- pd

replacePCA <- FALSE
if (replacePCA) {
	cat("*********\n")
	cat(" Replace genetic PCs made with merged-HM3 to those using just our samples\n")
	pd <- pd[,-which(colnames(pd) %in% paste("C",1:6,sep=""))]
	pData(MSet.genome) <- pd
	source("mergeWithGenetic.R")
	MSet.genome <- mergeWithGenetic(MSet.genome,snpPCFile)
	cat("Merge with new successful!\n")
	outF <- sprintf("%s/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",
		rootDir)
	save(MSet.genome, file=outF)
}

cat("* Running probe-wise analysis\n")
outDir <- sprintf("%s/dmp_DX_QTL/agesexPMI_PC12_%s",outRoot,dt)

dmpPref <- sprintf("%s/dmp..dmp_",outDir)
dmpFiles <- list(
	SEX=sprintf("%sSEXMale.topTable.%s.txt",dmpPref,dt),
	AGE=sprintf("%sAGE.topTable.%s.txt",dmpPref,dt),
	DX=sprintf("%sDXcase.topTable.%s.txt",dmpPref,dt)
)

if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)
outPref <- sprintf("%s/dmp.",outDir)

### Call DMPs
sink(sprintf("%s_%s.log",outPref,dt),split=TRUE)
tryCatch({
	print(table(pd$DX,pd$SEX))
	cat(sprintf("DMP finder: started at %s\n",format(Sys.time())))
	print(system.time(res <- dmpFinder_DX(MSet.genome,outPref,
						excludeTechReps=FALSE,incGeneticPCs=TRUE,
						incSlide=FALSE)))
	cat(sprintf("DMP finder: completed at %s\n",format(Sys.time())))
}, error=function(ex) { 
		print(ex)
},finally={
		sink(NULL)
})

#### ------------------------------------------------------
#### annotate DMPs
#### ------------------------------------------------------
dt <- format(Sys.Date(),"%y%m%d")

outDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/agesexPMI_PC12_171129"
dmpPref <- sprintf("%s/dmp..dmp_",outDir)
dmpFiles <- list(
	SEX=sprintf("%sSEXMale.topTable.%s.txt",dmpPref,dt),
	AGE=sprintf("%sAGE.topTable.%s.txt",dmpPref,dt),
	DX=sprintf("%sDXcase.topTable.171129.txt",dmpPref,dt)
)
logFile <- sprintf("%s_annotateDMP_%s",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
# for pathway analysis
genes <- read.delim(annoFiles[which(annoFiles[,1]=="Genes"),3],
	sep="\t",h=T,as.is=T)
gene_GR <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),name=genes[,5])

locs <- getLocations(MSet.genome)
options(scipen=10)
FDR_cutoff <- 0.3

for (nm in "DX") {
	cat(sprintf("%s\n",nm))
	dmp <- read.delim(dmpFiles[[nm]],sep="\t",h=T,as.is=T)
	colnames(dmp)[4] <- "pval"
	colnames(dmp)[5] <- "qval"

	# pathway analysis
	cat("* Pathway analysis\n")
	pThresh <- 0.05; 
	pathwayFile <- annoFiles[which(annoFiles[,1]=="Pathways"),3]
	logFile <- sprintf("%s/dmp%s_pathway_stats_p%1.3f_%s.txt",outDir,nm,
		 pThresh,format(Sys.time(),"%y%m%d.%H%M"))

	sink(logFile,split=TRUE)
	tryCatch({
		pathStats <- dmp_pathwayORA(MSet.genome,dmp,gene_GR,thresh=pThresh,
					pathwayFile=pathwayFile,outDir=outDir,
					gmtFDRcutoff=0.05,selMode="pvalue")
	},error=function(ex) {print(ex)
	},finally={
		sink(NULL)
	})


	samps <- pData(MSet.genome)$Sample_Name
	idx <- grep("-R[23456]",samps)
	noTechReps <- MSet.genome[,-idx] # remove tech duplicates

	# plot boxplots of betas
	cat("* Plotting boxplots of betas\n")
	pdf(sprintf("%s.DIST.DX.betas.Q0.05.pdf",
				sub(".txt","",dmpFiles[[nm]])),
		width=8,height=8)
	tryCatch({
		pd <- pData(noTechReps)
		tmp <- pd$DIST.DX
		tmp[grep("Control",tmp)] <- "CONT"
		tmp[grep("Schiz",tmp)] <- "SCZ"
		tmp[grep("Bipolar",tmp)] <- "BPD"
		pd$DX2 <- factor(tmp,levels=c("CONT","SCZ","BPD"))

		pData(noTechReps) <- pd
		dmp_plotProbes_groups(noTechReps,
			selProbes=rownames(dmp)[which(dmp$qval< 0.05)],
			group="DX2",baseline_group="CONT",group2="SEX",
			rmLegend=TRUE)
	},error=function(ex){
		print(ex)
	},finally={
		dev.off()
	})

	# write locations
	cat("* Writing probe locations\n")
	cur_loc <- annotateSig(MSet.genome, dmp, annoFiles,qCutoff=FDR_cutoff,
				outDir=outDir,fPrefix=nm)
	cur_loc <- as.data.frame(cur_loc)

	# plot pvalue of window around dmp
	dmp$NAME <- rownames(dmp)
	Qthresh <- 0.1
	pdf(sprintf("%s.localP_Q%1.2f.pdf",sub(".txt","",dmpFiles[[nm]]),
				Qthresh),width=11,height=8)
	tryCatch({
		dmp_plotLocalP(MSet.genome, dmp, 
					   rownames(dmp)[which(dmp$qval < Qthresh)],
					   winWidth=30*1000)
	},error=function(ex) {
			print(ex)
	},finally={
		dev.off()
	})

	### write IGF2-centered beta values
	
	dmp_plotLocalBeta(MSet.genome,"cg07096953",winWidth=2000,
				writeMelted=TRUE, outPref=sub(".txt","",dmpFiles[[nm]])) 
	pd <- pData(noTechReps)
	betas <- getBeta(noTechReps[c("cg07096953","cg02613624","cg26401390")])
	tmpFile=sprintf("%s/IGF2_probes_sampleDat.txt",dirname(dmpFiles[[nm]]))
	plotProbes_showSamp(betas,pd,tmpFile)

	cur_loc$NAME <- rownames(cur_loc)
	x <- merge(x=cur_loc,y=dmp,by="NAME")
	x <- x[order(x$qval),]

	outF <- sprintf("%s.FDR%1.2f.locations.txt", 
				sub(".txt","",dmpFiles[[nm]]),FDR_cutoff)
	write.table(x,file=outF,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}
},error=function(ex){
	print(ex)
},finally={
	dev.off()
})
	
