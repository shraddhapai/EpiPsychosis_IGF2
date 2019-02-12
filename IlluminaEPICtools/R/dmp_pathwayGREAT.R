#' prepare files for manual pathway analysis by GREAT
#'
#' @param mset (GenomicRatioSet) 
#' @param dmpIn (DataFrame) output of minfi::dmpFinder()
#' @param geneFile (char) path to bed file with gene definitions
#' @param geneDomain (numeric) max distance (kb) from TSS at which a probe is considered
#' within the gene's domain. Used for GREAT analysis. Should be large enough to include
#' any probe that would be part of a gene's domain by GREAT's calculation. Should be small
#' enough that the resulting bed files can be uploaded to the GREAT server.
#' (http://bejerano.stanford.edu/great/public/html/)
#' @param outDir (char) path to output dir
#' @export
dmp_pathwayGREAT <- function(mset,dmpIn,geneFile,cutoff=0.05,geneDomain=50,outDir) {
	dt <- format(Sys.Date(),"%y%m%d")
	
	bg <- rownames(dmpIn)
	fgpval <- rownames(dmpIn)[which(dmpIn$pval<cutoff)]
	# extract foreground, background, get locations, make to GR.
	loc <- getLocations(mset)
	loc_fg <- loc[which(names(loc)%in% fgpval)]
	loc_bg <- loc[which(names(loc)%in% bg)]
	cat(sprintf("%i foreground (p < cutoff), %i background\n", 
				length(fgpval),length(loc_bg)))
	
	options(scipen=10) # dont write genomic coords in sci notation.

	genes <- read.delim(geneFile,sep="\t",h=F,as.is=T)
	gene_GR <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),strand=genes[,6]) 
	gene_GR <- resize(gene_GR,fix="start",width=1)
	gene_GR <- resize(gene_GR,width=geneDomain*1000,fix="center")

	tmp <- subsetByOverlaps(loc_fg,gene_GR)
	cat(sprintf("\tGene domain = %1.1f kb\n", geneDomain))
	cat(sprintf("\tForeground: %i total > %i in gene domains\n", length(fgpval),length(tmp)))
	oF <- sprintf("%s/dmp_GREAT_fg_%s.txt",outDir,dt)
	tmp <- as.data.frame(tmp);
	write.table(tmp[,1:3],file=oF,sep="\t",col=F,row=F,quote=F)

	tmp <- subsetByOverlaps(loc_bg,gene_GR)
	cat(sprintf("\tBackground: %i total > %i in gene domains\n", length(loc_bg),length(tmp)))
	oF <- sprintf("%s/dmp_GREAT_bg_%s.txt",outDir,dt)
	tmp <- as.data.frame(tmp);
	write.table(tmp[,1:3],file=oF,sep="\t",col=F,row=F,quote=F)

}
