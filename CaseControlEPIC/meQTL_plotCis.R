# plot summary stats on meQTL
rm(list=ls())
require(GenomicRanges)

resFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/meQTL/180110/regular_lmer_genoDisease0/withC1C2/eQTL.results.180110.txt"
rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files"
snpFile <- sprintf("%s/CaseControlSNP/meQTL_postimpute/171207/geno.meQTL.snps.171207.traw",rootDir)
mDataFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata"
# genotypes

outDir <- dirname(resFile)

load(mDataFile)
res <- read.delim(resFile,sep="\t",h=T,as.is=T)
MSet.genome <- MSet.genome[unique(res$cpg)]
mLocs <- getLocations(MSet.genome)

cpos <- regexpr(":",res$snps)
chr <- paste("chr",substr(res$snps,1,cpos-1),sep="")
pos <- as.integer(substr(res$snps,cpos+1,nchar(res$snps)))
snpGR <- GRanges(chr,IRanges(pos,pos))
names(snpGR) <- res$snps

mpos <- start(mLocs[res$cpg])
spos <- start(snpGR[res$snps])
dsnp <- mpos-spos

res$dsnp <- dsnp/1e6
dpos <- regexpr(":",res$snps)
res$chrom <- as.integer(substr(res$snps,1,dpos-1))

res <- subset(res,FDR < 0.05)
outDir <- dirname(resFile)
outFile <- sprintf("%s/snpd_vs_meth.pdf", outDir)
pdf(outFile)
par(las=1,bty='n',cex.axis=1.3)
plot(res$dsnp, abs(res$slope*100), col=rgb(0,0,0,0.5), cex=0.5,pch=16,
	xlab="SNP-CpG distance (Mb)", 
	ylab="% methylation change\nper minor allele")
abline(v=0,lty=3,col='red')
# letters for chroms
plot(res$dsnp, (res$slope*100), col=rgb(0,0,0,0.5), cex=0.5,pch=16,
	xlab="SNP-CpG distance (Mb)", 
	ylab="% methylation change\nper minor allele",type='n')
letters <- strsplit("abcdefghijklmnopqrstuvw","")[[1]]
text(res$dsnp, (res$slope*100), labels=letters[res$chrom],cex=0.4)
abline(v=0,lty=3,col='red')
dev.off()

