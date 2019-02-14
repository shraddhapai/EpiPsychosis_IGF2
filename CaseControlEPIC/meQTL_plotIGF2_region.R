rm(list=ls())
require(GenomicRanges)
require(Gviz)
require(GenomicInteractions)

rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files"
resFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/meQTL/180110/regular_lmer_genoDisease0/withC1C2/allOutput.Rdata"
load(resFile)
mData <- sprintf("%s/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",
				 rootDir)
load(mData)
require(minfi)
mLocs <- getLocations(MSet.genome)

snpFile <- sprintf("%s/CaseControlSNP/meQTL_postimpute/171207/geno.meQTL.snps.171207.traw",rootDir)
# which type of snp is it? pgc associated/credible/near dmp
source("meQTL_prepareGeno.R")
snpDat <- meQTL_prepareGeno(snpFile=snpFile)
snpGR <- snpDat$snpGR
snps <- snpDat$geno

cat(sprintf("Got %i meQTL associations\n", nrow(out)))
bonfp <- 0.05/nrow(out)

tmp <- out[which(out$FDR<0.05),]
tmp$PASS_BONFERRONI <- tmp$p < bonfp
outF <- sprintf("%s/FDRpass.txt", dirname(resFile))
write.table(tmp,file=outF,sep="\t",col=T,row=F,quote=F)

dmpGR <- GRanges("chr11",IRanges(2154255,2154952))
dmpGR <- resize(dmpGR, fix="center",width=1e6)
idx <- findOverlaps(snpGR, dmpGR)
mysnps <- snpGR$name[queryHits(idx)]

out2 <- subset(out, snps %in% mysnps)
cat(sprintf("%i snps in IGF2 area\n", nrow(out2)))
out2$snps  <- as.character(out2$snps)
outDir <- dirname(resFile)

pdf(sprintf("%s/IGF2_meQTL.pdf",outDir),width=11,height=5)
for (k in unique(out2$cpg)) {
print(k)
	cur  <- subset(out2, cpg %in% k)
	mypos <- start(mLocs[k])
	loc <- as.integer(unlist(lapply(strsplit(cur$snps,":"),function(x) x[2])))
	cur$loc <- loc
	cur$d <- cur$loc-mypos
	cur <- cur[order(cur$d),]
	plot(cur$d/1000,-log10(cur$p),type='o',pch=16,main=k,
		xlab="distance from probe (Kb)",cex=0.7)
	idx <- which(cur$FDR < 0.05)
	if (any(idx)) {
		points(cur$d[idx]/1000,-log10(cur$p[idx]),col='red',pch=16,cex=0.7)
	}
	# closeup.
	plot(cur$d/1000,-log10(cur$p),type='o',pch=16,main=k,
		xlab="distance from probe (Kb)",cex=0.7,xlim=c(-100,100))
	idx <- which(cur$FDR < 0.05)
	if (any(idx)) {
		points(cur$d[idx]/1000,-log10(cur$p[idx]),col='red',pch=16,cex=0.7)
	}
}
dev.off()

source("outdated/plotEQTL.R")
out2$cpg <- as.character(out2$cpg)
mvals <- getBeta(MSet.genome)
pheno <- pData(MSet.genome)
for (k in unique(out2$cpg)) {
	print(k)
	idx <- which(out2$cpg %in% k & out2$FDR < 0.05)
	if (any(idx)) {
		pdf(sprintf("%s/%s_genoMeth.pdf",outDir,k),width=11,height=8)
		suppressWarnings(plotEQTL(out2[idx,,drop=FALSE],
			snps,mvals,pheno,snpLocs))
		dev.off()
	}
}


