#' utility to standardize geno prep for meQTL
#' @param snpFile (char) .traw file
#' @return (list) (1) snpGR: snp GRanges, (2) snpinfo: first 6 cols
#' (3) geno : genotype matrix
meQTL_prepareGeno <- function(snpFile) {
	snps <- read.delim(snpFile,sep="\t",h=T,as.is=T)
	sampNames <- sub("X","",colnames(snps))[-(1:6)]
	upos <- regexpr("_",sampNames)
	sampNames <- substr(sampNames,upos+1,nchar(sampNames))
	sampNames <- sub("\\.","-",sampNames)

	snpinfo <- snps[,1:6]
	snpGR <- GRanges(paste("chr",snpinfo$CHR,sep=""),
		IRanges(snpinfo$POS,snpinfo$POS+1),name=snpinfo$SNP)
	rm(snpinfo)
	snps <- snps[,-(1:6)]
	colnames(snps) <- sampNames; rm(sampNames)
	rownames(snps) <- snpGR$name
	toomany_na <- which(rowSums(is.na(snps))>5)
	if (any(toomany_na)) {
		cat(sprintf("These snps have >5 NAs and are being excluded\n"))
		cat(sprintf("%s\n",paste(snpGR$name[toomany_na],collapse=",")))
		snps <- snps[-toomany_na,]
	}
# we know 90 and 95 are technical replicates with duplicate identifiers
# assigned by the tissue bank.
colnames(snps)[which(colnames(snps)=="90-R1")] <- "95-R4"
colnames(snps)[which(colnames(snps)=="90-R2")] <- "95-R5"
colnames(snps)[which(colnames(snps)=="90-R3")] <- "95-R6"

out <- list(snpGR=snpGR,geno=snps)
return(out)

}
