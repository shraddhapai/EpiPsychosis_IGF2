#' get granges of genes
#'
#' @param gFile (char) path to gene definition file
#' @return (GRanges)
#' @import GenomicRanges
#' @export
getGenes <- function(gFile) {
	dat <- read.delim(gFile,sep="\t",h=T,as.is=T)
	cList <- c(paste("chr",1:22,sep=""),"chrX","chrY")
	dat <- subset(dat,chrom %in% cList)

	gene_GR <- GRanges(dat$chrom,IRanges(dat$txStart,dat$txEnd),
		strand=dat$strand)
	gene_GR$name	<- dat$name
	gene_GR$name2	<- dat$name2

	gene_GR
}
