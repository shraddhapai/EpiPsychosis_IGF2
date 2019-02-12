rm(list=ls())


matchGeno <- function(rnaFile,genoFile,outFile) {
	rna <- read.table(rnaFile,sep="\t",h=T,as.is=T)
	baseF <- basename(rnaFile)
	upos <- gregexpr("_", baseF)[[1]]
	baseF <- substr(baseF,upos[1]+1,upos[2]-1)
	baseF <- as.integer(baseF)
	cat(sprintf("Sample %s\n", baseF))
	
	geno <- read.table(genoFile,sep="\t",h=T,as.is=T)
	genoinfo <- geno[,1:6]; geno <- geno[,-(1:6)]
	
	sampName <- colnames(geno); 
	sampName <- sub("X","",sampName)
	upos <- regexpr("_", sampName)
	sampName <- substr(sampName,upos+1,nchar(sampName))
	idx <- grep(".R[23]",sampName)
	geno <- geno[,-idx]; sampName <- sampName[-idx]
	sampName <- as.integer(sub(".R1","",sampName))
	
	# match snps
	genosnp <- paste(genoinfo$CHR,genoinfo$POS,sep="_")
	idx <- which(!duplicated(genosnp))
	geno <- geno[idx,]; genoinfo <- genoinfo[idx,]; genosnp <- genosnp[idx]
	
	rnasnp <- paste(rna$CHR,rna$POS,sep="_")
	common <- intersect(genosnp,rnasnp)
	
	idx <- which(genosnp %in% common)
	geno <- geno[idx,];genoinfo <- genoinfo[idx,];genosnp <- genosnp[idx];rm(idx)
	idx <- which(rnasnp %in% common)
	rna <- rna[idx,]; rnasnp <- rnasnp[idx]; rm(idx)
	cat(sprintf("%i geno and %i rna left \n",nrow(geno),nrow(rna)))
	
	
	midx <- match(genosnp, rnasnp)
	if (all.equal(rnasnp[midx],genosnp)!=TRUE) {
		cat("rna-geno not match"); browser()
	}
	rna <- rna[midx,]
	
	# same allele must be major/minor
	idx <-which(genoinfo$COUNTED == rna$COUNTED & genoinfo$ALT == rna$ALT)
	geno <- geno[idx,]; rna <- rna[idx,]
	cat(sprintf("%i after matching genotype to rna\n",nrow(geno)))
	rm(midx,idx,genoinfo,genosnp,rnasnp)
	
	# look up geno
	gidx <- which(sampName==baseF)
	if (is.na(gidx)) {
		cat(sprintf("*** SKIPPING sample %s, match not found!\n", baseF))
	} else {
	# count allele match
	g <- geno[,gidx]
	naidx <- which(is.na(g) | is.na(rna[,7]))
	if (any(naidx)) { g <- g[-naidx]; rna <- rna[-naidx,] }
	tot <- sum(g==rna[,7])
	cat(sprintf("%s: %i of %i alleles match (%1.2f)\n", 
		baseF,tot, nrow(geno),(tot/nrow(rna))))
	cat(sprintf("%s\t%i\t%i\t%1.2f%%\n",baseF,tot,nrow(geno),(tot/nrow(rna))),
		file=outFile,append=TRUE)
	}
	
}

genoFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/final_hg38/SUS19399.hg19.sorted_-clean.hg38.rnamatch.traw"
rnaDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlRNAseq/var_calls"

dt <- format(Sys.Date(),"%y%m%d")
logFile <- sprintf("matchGeno_%s.log",dt)
sink(logFile,split=TRUE)

tryCatch({
	fNames <- dir(rnaDir,pattern="vcf.plink.traw")
	cat(sprintf("Got %i files\n",length(fNames)))
	outFile <-"rna2snp_matches.txt"
	system(sprintf("cat /dev/null > %s\n",outFile))
	cat("Sample\tsnps_matched\tsnps_tested\tpct\n",file=outFile,append=TRUE)	

	for (curF in fNames) {
		rnaFile <- sprintf("%s/%s",rnaDir,curF)
		cat("-----\n")
		matchGeno(rnaFile=rnaFile,genoFile=genoFile,outFile)	
		cat("-----\n")
	}
	
	dat <- read.delim("rna2snp_matches.txt",sep="\t",h=T,as.is=T)
	pdf("rna2snp_density.pdf")
	dat[,4] <- as.numeric(sub("%","",dat[,4]))#
	plot(density(dat[,4]), cex.axis=1.5,las=1,bty='n')
	dev.off()
	
},error=function(ex){
	print(ex)
},finally={
	cat("closing")
	sink(NULL)
})



