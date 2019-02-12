#' meQTL for IGF2-TH interactions, specifically testing for geno*dx 
rm(list=ls())
rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files"
snpFile <- sprintf("%s/CaseControlSNP/meQTL_postimpute/171207/geno.meQTL.snps.pruned.171207.traw",rootDir)
dmpFile <- sprintf("%s/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/agesexPMI_PC12_171129_genes.txt",
				   rootDir)
mData <- sprintf("%s/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",
				 rootDir)

outRoot <- sprintf("%s/meQTL",dirname(dmpFile))
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s/IGF2TH_genoDisease1",outRoot,dt)
if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)

logFile <- sprintf("%s/IGF2_genoDx_%s.log", outDir, dt)
sink(logFile,split=TRUE)
tryCatch({
cpg2test <- c("cg02613624","cg22956483","cg07096953")
#lim2snp <- c("11:2186335","11:2190519")

load(mData)
MSet.genome <- MSet.genome[cpg2test]
mLocs <- getLocations(MSet.genome)

mvals <- getBeta(MSet.genome)
# average probes in DMR2 region
mvals <- t(as.matrix(colMeans(mvals)))
rownames(mvals) <- "IGF2_DMR"
tmp <- GRanges("chr11",IRanges(min(start(mLocs)),max(end(mLocs))))
mLocs <- tmp
names(mLocs) <- "IGF2_DMR"

pheno <- pData(MSet.genome)
colnames(mvals) <- pheno$Sample_Name

source("eQTL_lmer.R")
out <- eQTL_lmer(snpFile,mvals,mLocs,pheno,
	addGenoDiseaseTerm=TRUE,MAX_DIST=1e6)
write.table(out$res,file=sprintf("%s/pruned_pvalue.stats.txt",outDir),sep="\t",col=T,
	row=F,quote=F)

res <- out$res
cat(sprintf("Num tests = %i ; Num SNPs= %i; Num CpGs=%i\n",
	nrow(res),length(unique(res$snp)), length(unique(res$cpg))))
cat("FDR summary\n")
print(summary(res$FDR))

},error=function(ex){print(ex)},finally={sink(NULL)})

cpos <- regexpr(":", res$snps)
pos <- as.integer(substr(res$snps,cpos+1,nchar(res$snps)))

dsnp <- pos - start(mLocs)[1]




