#' plot locus-wise coverage across samples
#' 
#' 1. Plots base-level coverage. Strand-split. Strands pooled.
#' 2. Plots target-level coverage. 
rm(list=ls())
source("multiplot.R")
source("poolStrands.R")
source("parseGTF.R")
require(reshape2)
require(ggplot2)
require(GenomicRanges)
inDir <- "/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/output/methylation"
geneFile <- "/home/shraddhapai/genome_annotation/grch38/gencode.v26.annotation.gtf.genesonly.txt"

target_groups <- c("Adol_Reg","NEURODEV","scz2_plt0.0001.CTCF",
	"scz2_plt0.0001.EnhProm","AdolEPIC") #,"chrX")

for (gp in target_groups[2:4])  {
	fList <- dir(inDir, pattern=sprintf("%s.records.Rdata",gp))
	isEmpty <- matrix(NA,nrow=length(fList),ncol=2)
	ctr <- 1
	isEmptyName <- list() # which loci have no coverage
	for (fName in fList) {
		sampName <- sub(sprintf(".%s.records.Rdata",gp),"",fName)
		out <- poolStrands(sprintf("%s/%s",inDir,fName),getTargetGR=TRUE)
		rec <- out$rec
		names(rec) <- names(out$target_GR)

		empty <- c()
		for (k in 1:length(rec)) {
			if (is.null(dim(rec[[k]]))) empty <- c(empty,k)
		}
		isEmpty[ctr,] <- c(length(empty),length(rec))
		isEmptyName[[sampName]] <- out$target_GR[empty]
		ctr <- ctr+1
	}
	
	empty_name <- unique(unlist(lapply(isEmptyName, function(x) x$name)))
	empty_mat <- matrix(0,nrow=length(empty_name),ncol=length(fList))
	rownames(empty_mat) <- empty_name
	for (k in 1:length(isEmptyName)) {
		cur_empty<-isEmptyName[[k]]$name
		empty_mat[which(rownames(empty_mat)%in% cur_empty),k]<-1
	}

	target_GR <- out$target_GR
	tot <- rowSums(empty_mat)
	idx <- which(tot == ncol(empty_mat))
	cat(sprintf("%i targets have no coverage in all samples\n",length(idx)))
	tmp <- target_GR[which(target_GR$name %in% rownames(empty_mat)[idx])]	
	df <- as.data.frame(tmp)
	df$pctEmpty <- 100

	idx <- which(tot > round(ncol(empty_mat)/2))
	cat(sprintf("%i targets have no coverage in half samples\n",length(idx)))
	tmp <- target_GR[which(target_GR$name %in% rownames(empty_mat)[idx])]	
	df2 <- as.data.frame(tmp)	
	df2 <- df2[-which(df2$name %in% df$name),] # remove duplicates
	df2$pctEmpty <- 50
	
	df <- rbind(df,df2)

	# find the nearest genes to the sites with missing coverage
	cat("Mapping missing targets to nearest genes\n")
	genes <- read.delim(geneFile,sep="\t",h=F,as.is=T)
	gene_GR <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),name=genes[,6])
	colnames(genes) <- paste("gene",colnames(genes),sep=".")
	df_GR <- GRanges(df$seqnames,IRanges(df$start,df$end),name=df$name)
	d <- distanceToNearest(df_GR, gene_GR)
	df2gene <- cbind(df[queryHits(d),], genes[subjectHits(d),])
	write.table(df2gene,file=sprintf("%s.noCvg.txt",gp),sep="\t",col=T,row=F,quote=F)

	rownames(isEmpty) <-fList
	colnames(isEmpty) <- c("num_targets_nocvg","total_targets")
	write.table(isEmpty,file=sprintf("%s.isEmpty.txt",gp),sep="\t",col=T,row=T,
		quote=F)
}
