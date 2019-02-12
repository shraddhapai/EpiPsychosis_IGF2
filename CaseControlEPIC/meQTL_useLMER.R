#' meQTL for SNPs around Case-control DMPs
#' use LMER
rm(list=ls())
require(minfi)
require(GenomicRanges)

# ------------------
# input files
rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files"

#### genetic
# genotypes
snpFile <- sprintf("%s/CaseControlSNP/meQTL_postimpute/171207/geno.meQTL.snps.171207.traw",rootDir)
# which type of snp is it? pgc associated/credible/near dmp
snpSrc <- sprintf("%s/CaseControlSNP/meQTL_postimpute/171207/geno.meQTL.snpinfo.171207",rootDir)
# location of snps
snpLocFile <- sprintf("%s/CaseControlSNP/plinkQC/postimpute/postimputeqc.bim",
					  rootDir)
#### methylation
mData <- sprintf("%s/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",
				 rootDir)
dmpFile <- sprintf("%s/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/agesexPMI_PC12_171129_genes.txt",
				   rootDir)
plinkMds <- sprintf("%s/CaseControlSNP/SUS19399_New/hg19/SUS19399.hg19.sorted_-clean_PCA.eigenvec",
	rootDir)

### merge with nature of SNP . credible SNP2, or index SNPs or in dmp zone.
currModeTypes <- "excludeCisDMP" #c("regular") # regular for cis ; excludeCisDMP for trans
FDRthresh <- 0.05
outRoot <- sprintf("%s/meQTL",dirname(dmpFile))
# ------------------------------------------------
dt <- format(Sys.Date(),"%y%m%d")

for (modeType in currModeTypes) {
	for (addGenoDiseaseTerm in FALSE) { #c(FALSE,TRUE)) {
		outDir <- sprintf("%s/%s/%s_lmer_genoDisease%i",outRoot,dt,modeType,
						  as.integer(addGenoDiseaseTerm))
		if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)
		
		logFile <- sprintf("%s/meQTL_dmp_%s.log",outDir,dt)
		sink(logFile,split=TRUE)
		tryCatch({
		
		cat("* Read DMP locs\n")
		dmp <- read.delim(dmpFile,sep="\t",h=TRUE,as.is=T)
		dmp <- subset(dmp, z_sidak_p < FDRthresh)
		cat(sprintf("\t%i DMPs with Q < %1.2f\n", nrow(dmp),FDRthresh))
		dmp_GR <- GRanges(dmp$seqnames,IRanges(dmp$start,dmp$end),
						  name=dmp$gene)
		cat(sprintf("\tGot %i DMP regions\n",length(dmp_GR)))
		
		cat("* Read methylation data\n")
		load(mData)
		locs <- getLocations(MSet.genome)
		ol <- findOverlaps(locs,dmp_GR) # limit to CpG in DMP regions
		idx <- unique(queryHits(ol))
		cat(sprintf("\t%i probes in DMP regions\n",length(idx)))
		mvals <- getBeta(MSet.genome[idx,])
		mLocs <- getLocations(MSet.genome[idx,])
		pheno <- pData(MSet.genome)
		colnames(mvals) <- pheno$Sample_Name

		cat("* Read SNP source\n")
		info <- read.delim(snpSrc,sep="\t",h=F,as.is=T)
		colnames(info) <- c("snps","snp.source")
		info <- info[!duplicated(info),]
		if (modeType %in% "excludeCisDMP") {
			cat("*** Limiting to PGC2 snps:")
			idx <- grep("SNPs in DMP", info$snp.source)
			lim2snp <- info$snps[-idx]
			cat(sprintf("%i SNPs\n", length(lim2snp)))
		} else {
			lim2snp <- "*"
		}

		
		# ------------------------------------------------
		cat("* Run MEQTL -- lm\n")
		source("eQTL_lmer.R"); 
		pdf(sprintf("%s/meQTL_pvalue_hist_%s.pdf", outDir,dt))
		tryCatch({
			# if testing cis-eQTL impose distance constraint
			if (modeType=="regular") MAX_DIST <- 1e6 else MAX_DIST <- Inf
			out <- eQTL_lmer(snpFile,mvals,mLocs,pheno,lim2snps=lim2snp,
						   addGenoDiseaseTerm=addGenoDiseaseTerm,
								MAX_DIST=MAX_DIST)
		}, error=function(ex) {
			print(ex)
		},finally={
			dev.off()
		})
		rm(mvals,pheno)
		
		snps 	<- out$snps
		mvals 	<- out$M
		pheno	<- out$cvrt
		out		<- out$res
		
		# process results
		outFile <- sprintf("%s/eQTL.results.%s.txt",outDir,dt)
		options(scipen=10)
		out <- out[order(out$FDR),]
		write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)
		res <- out
		colnames(res)[1] <- "gene"
		pdf(sprintf("%s/phist.pdf",outDir))
		hist(out$p,n=100,main="meQTL pvalues"); dev.off()
		save(out,file=sprintf("%s/allOutput.Rdata",outDir))

		if (sum(res$FDR < 0.05)>0) {
			res <- subset(res, FDR < 0.05)
			write.table(res,file=sprintf("%s/passFDR.txt", outDir),
				sep="\t",col=T,row=F,quote=F)
		}
		
		cat("Process completed successfully\n")
		},error=function(ex) {
			print(ex)
		},finally={
			cat("Closing log\n")
			sink(NULL)
		})
  }
}
