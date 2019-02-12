#' FET for various neurobiological genesets 
#' 
#' @param diffEx (data.frame) diffex output from edgeR glmLRTc
diffEx_ORA <- function(diffEx) {

source("../../util/readPathwayFile.R")
source("../../util/cleanPathwayName.R")

# BrainSpan cell dev phases
braindev <- "/Users/shraddhapai/Google Drive/genome_annotation/BrainSpan/BrainSpan_TableS13.gmt"
cellGeneDir <- "/Users/shraddhapai/Google Drive/genome_annotation/BrainSpan/Darmianis_TableS3"
mir483_targets <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/HanZoghbi_GenesDev_2013_TableS2.txt"
lakeFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/LakeZhang_2016_Science/Lake_TableS5_in.txt"

# -------------------------------------
# compile pathway list
braindev <- readPathways(braindev,MIN=0,MAX=1000)
# Darmanis Table S3 - celltype specific
fList <- dir(cellGeneDir,pattern="txt")[-6]
celltypes <- list()
cdir <- cellGeneDir
for (fName in fList) {	
   dpos <- regexpr("\\.",fName)
   curtype <- substr(fName,dpos+1,nchar(fName)-4)
   dat <- read.delim(sprintf("%s/%s",cdir,fName),h=F,
                     as.is=TRUE)[,1]
   celltypes[[curtype]] <- dat
   #cat(sprintf("\t%s:%i markers\n",curtype,length(dat)))
   
}

# neuronal subtypes from Lake et al.
neurogene <- read.delim(lakeFile,sep="\t",h=T,as.is=T)
neuroset <- list()
for (k in unique(neurogene$cluster)) {
	curgroup <- neurogene[which(neurogene$cluster==k),2]
	if (length(curgroup)>=10) {
		neuroset[[k]] <- curgroup
	} 
}

dat <- read.delim(mir483_targets,sep=" ",h=T,as.is=T)
mir483 <- list(mir483_targets=dat[,1])

allpath <- c(braindev,celltypes,neuroset,mir483)

# --------------
# get diffEx genes
is_sig <- diffEx$gene.name[which(diffEx$PValue < 0.05)]
not_sig <- setdiff(diffEx$gene.name,is_sig)

res <- matrix(NA,nrow=length(allpath), ncol=6)
ctr <- 1
fg_genes <- c()
for (nm in names(allpath)) {
	curgroup <- allpath[[nm]]
	fg_in_group <- intersect(is_sig, curgroup)
	bg_in_group <- intersect(not_sig,curgroup)
	mat <- matrix(NA,nrow=2,ncol=2)
	mat[1,1] <- length(fg_in_group)
	mat[1,2] <- length(is_sig)-length(fg_in_group)
	mat[2,1] <- length(bg_in_group)
	mat[2,2] <- length(not_sig)- length(bg_in_group)
	x <- fisher.test(mat)
	#cat(sprintf("%s: %i total: %i in sig (p< %1.2e)\n",
	#		nm,length(curgroup),length(fg_in_group),x$p.value))
	
	res[ctr,] <- c(length(curgroup),length(fg_in_group),length(bg_in_group),
		length(fg_in_group)/length(is_sig), length(bg_in_group)/length(not_sig),
		x$p.value)
	fg_genes <- c(fg_genes,paste(fg_in_group,collapse=","))
	
	ctr <- ctr+1
}
colnames(res) <- c("set_size","FG overlap","BG overlap",
		"%FG","%BG","FET.pvalue")
res[,4:5] <- round(res[,4:5]*100,digits=2)
res <- as.data.frame(res)
res$FG_genes <- fg_genes
res$Geneset_name <- names(allpath)
res <- res[,c(ncol(res),1:(ncol(res)-1))]
res <- res[order(res$FET.pvalue),]
cat(sprintf("Bonferroni p is: %1.2e\n", 0.05/nrow(res)))

return(res)
}

