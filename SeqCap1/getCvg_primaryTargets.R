#' plot coverage of primary targets
rm(list=ls())

primaryFile <- "/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/input/OID43737_hg19_161128_design_deliverables/OID43737_hg19_161128_primary_targets.hg38.bed"
inDir <- "/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/output/align/newPicard"
outDir <- "/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/output/methylation"

# -----------------------------------------------
require(GenomicRanges)
source("getM_GRanges.R")

if (!file.exists(outDir)) dir.create(outDir)

# Label target groups
target_groups <- "CaseControl" #c("Adol_Reg","NEURODEV","scz2_plt0.0001.CTCF",
	#"scz2_plt0.0001.EnhProm","AdolEPIC","CaseControl") #,"chrX")
targets <- read.delim(primaryFile,sep="\t",h=F,as.is=T)
mygroup <- rep("",nrow(targets))
for (g in target_groups) {
	mygroup[grep(g,targets[,4])] <- g
}
targets$group <- mygroup
targets <- subset(targets, group %in% c("scz2_plt0.0001.CTCF",
	"scz2_plt0.0001.EnhProm","NEURODEV","CaseControl"))
targets <- subset(targets, targets[,1] %in% paste("chr",c(1:22,c("X","Y","M")),
	sep=""))
dt <- format(Sys.Date(),"%y%m%d")

# --------------------------------------
# begin work
logFile <- sprintf("%s/getCoverage_%s.log", outDir,dt)
sink(logFile,split=TRUE)
tryCatch({

fList <- dir(path=inDir,pattern="methylation_results.txt.gz$")
cat(sprintf("Got %i samples\n",length(fList)))

for (samp in fList) {  #loop over samples
	inFile <- sprintf("%s/%s",inDir,samp)
	cat("---------------------------\n")
	cat(sprintf("%s\n",samp))
	cat("---------------------------\n")
	baseF <- sub(".methylation_results.txt.gz","",basename(samp))
	for (g in unique(targets$group)) { # loop over types of target regions
		cat(sprintf("\t%s\n",g))
		currt <- subset(targets, group %in% g)
		gr <- GRanges(currt[,1],IRanges(currt[,2],currt[,3]),name=currt[,4])
	
		# get methylation counts related to each target
		rec <- getRec_GRanges(inFile,gr,verbose=FALSE,numCores=16)
		names(rec) <- currt[,4]

		# high coverage
		# get methylation counts related to each target
		y <- lapply(rec, function(x) {
			if (length(dim(x))<2) return(cbind(C=0,COV=0,pctM=NA))
			x <- subset(x, CT_count>=10) # "high confidence"
			if (length(dim(x))<2) return(cbind(C=0,COV=0,pctM=NA))
		    
			#if (is.na(x)) return(cbind(C=0,COV=0,pctM=NA))
			#else {
			cbind(C=sum(x$C_count),COV=sum(x$CT_count),
				pctM=1-(sum(x$C_count)/sum(x$CT_count)))}
	)
			#})
		names(y) <- gr$name
		z2 <- data.frame(do.call("rbind",y))
		z2$pctM <- round(z2$pctM,digits=3)
		z2$target_length <- width(gr)
		z2$name <- names(y)
		rm(y)
	
		outFile <- sprintf("%s/%s.%s.records.Rdata", outDir,baseF,g)
		target_gr <- gr
		save(rec,target_gr,file=outFile)

#		write.table(z,file=sprintf("%s/%s.%s.summary.txt",outDir,baseF,g),
#			sep="\t",col=TRUE,row=F,quote=F)
		write.table(z2,file=sprintf("%s/%s.%s.cvgGe10.summary.txt",outDir,baseF,g),
			sep="\t",col=TRUE,row=F,quote=F)
		cat("\n")
	}
}
},error=function(ex) {
	print(ex)
}, finally={
	sink(NULL)
})
