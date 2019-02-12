# plot methylation with changing genotype
rm(list=ls())
require(GenomicRanges)

rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files"
snpFile <- sprintf("%s/CaseControlSNP/meQTL_postimpute/171207/geno.meQTL.snps.171207.traw",rootDir)
mData <- sprintf("%s/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",rootDir)
resFile <- sprintf("%s/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/meQTL/180110/regular_lmer_genoDisease0/withC1C2/eQTL.results.180110.txt",rootDir)

probeList <- c("cg07096953","cg02613624","cg22956483","cg26401390")

res <- read.delim(resFile,sep="\t",h=T,as.is=T)
res <- subset(res,FDR <= 0.1)
res <- subset(res,cpg %in% probeList)

source("meQTL_prepareGeno.R")
snpDat <- meQTL_prepareGeno(snpFile=snpFile)
gr <- snpDat$snpGR
geno <- snpDat$geno

idx <- which(gr$name %in% res$snps)
gr <- gr[idx]
geno <- geno[idx,]

# remove tech reps
load(mData)
MSet.genome <- MSet.genome[unique(res$cpg)]
pd <- pData(MSet.genome); M <- getBeta(MSet.genome)
idx <- grep("-R[23456]",pd$Sample_Name)
pd <- pd[-idx,]; M <- M[,-idx]
colnames(M) <- pd$Sample_Name

idx <- grep("-R[23456]",colnames(geno))
geno <- geno[,-idx]

common <- intersect(pd$Sample_Name, colnames(geno))
geno <- geno[,which(colnames(geno) %in% common)]
pd <- pd[which(pd$Sample_Name %in% common),]
M <- M[,which(colnames(M) %in% common)]

pd <- pd[order(pd$Sample_Name),]
geno <- geno[,order(colnames(geno))]
M <- M[,order(colnames(M))]

if (all.equal(colnames(M),colnames(geno))!=TRUE) {cat("M g don't mathc"); browser()}
if (all.equal(colnames(M),pd$Sample_Name)!=TRUE) {cat("M pd don't match"); browser()}

pd <- as.data.frame(pd)
pd$DX <- factor(pd$DX,levels=c("case","control"))
plotList <- list()
for (k in 1:nrow(res)) {
	pd_tmp <- pd
	pd_tmp$geno <- factor(as.integer(geno[res$snps[k],]))
	pd_tmp$M <- M[res$cpg[k],]
	p <- ggplot(pd_tmp,aes(geno,M))+geom_boxplot(aes(colour=DX))
	p <- p + geom_point(position=position_dodge(width=0.75),aes(group=DX),
		cex=0.5)
	p <- p + ggtitle(sprintf("%s:%s, Q < %1.2f", res$snps[k], res$cpg[k],
		res$FDR[k]))
	plotList[[k]] <- p
}
source("multiplot.R")
outDir <- dirname(resFile)
pdfFile <- sprintf("%s/cismeQTL_IGF2.pdf",outDir)
pdf(pdfFile)
multiplot(plotlist=plotList, layout=matrix(1:4,ncol=2))
dev.off()

