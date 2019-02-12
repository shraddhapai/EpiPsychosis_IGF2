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

###out <- out[which(out$p < bonfp),]
###cat(sprintf("\t%i associations pass Bonferroni p (p < %1.2e)\n",
###	nrow(out),bonfp))
###outF <- sprintf("%s/bonferroniPass.txt", dirname(resFile))
###write.table(out,file=outF,sep="\t",col=T,row=F,quote=F)
# location of snps


dmpGR <- GRanges("chr11",IRanges(2154255,2154952))
dmpGR <- resize(dmpGR, fix="center",width=1e6)
idx <- findOverlaps(snpGR, dmpGR)
mysnps <- snpGR$name[queryHits(idx)]

out2 <- subset(out, snps %in% mysnps)
cat(sprintf("%i snps in IGF2 area\n", nrow(out2)))
out2$snps  <- as.character(out2$snps)
outDir <- dirname(resFile)

#### ---------------------------------------
#### Gviz plot of interactions
###trackSaveDir <- "/home/shraddhapai/Epigenetics/NARSAD/geneModels"
### 
#### assemble granges for plot
###out2 <- out2[which(out2$FDR < 0.05),]
###mLocDF <- as.data.frame(mLocs[unique(out2$cpg)])
###mLocDF <- mLocDF[,c("seqnames","start","end")]
###out2$cpg <- as.character(out2$cpg)
###midx <- match(out2$cpg, rownames(mLocDF))
###if (all.equal(rownames(mLocDF)[midx],out2$cpg)!=TRUE) {
###	cat("cpg names don't match\n"); browser()
###}
###out2$cpg_chrom <- mLocDF$seqnames[midx]
###out2$cpg_pos <- mLocDF$start[midx]
###
###snpLoc <- out2$snps; tmp <- strsplit(snpLoc,":")
###tmp_chrom <- unlist(lapply(tmp,function(x)x[1]))
###tmp_pos <- unlist(lapply(tmp,function(x)x[2]))
###snpLoc <- data.frame(chr=tmp_chrom,pos=as.integer(tmp_pos))
###
###ylim <- c(min(c(mLocDF$start,snpLoc$pos)), max(c(mLocDF$start, snpLoc$pos)))
###
###tSaveFile <- sprintf("%s/IGF2_chr11_2000000_2300000.Rdata",trackSaveDir)
###cat("Loading from file!\n")
###load(tSaveFile)
###displayPars(refGenes) <- list(showId=TRUE,fill="darkgreen",
###			  cex=0.5,col=NULL)
#### prepare for plot view
###gr1 <- GRanges(out2$cpg_chrom,IRanges(out2$cpg_pos,out2$cpg_pos))
###gr2 <- GRanges(paste("chr",snpLoc$chr, sep=""),IRanges(snpLoc$pos,snpLoc$pos))
###ints		<- GenomicInteractions(gr1,gr2,counts=25)
###intTrack	<- InteractionTrack(name="meQTL",ints,
###							chromosome="chr11") 
###displayPars(intTrack) <- list(plot.trans=FALSE, 
###								  plot.outside=TRUE,
###								  col.outside="lightblue",
###			interaction.measure=out2$FDR,
###			interaction.dimension.transform="log")
###
###tmp <- out2[!duplicated(out2$cpg_pos),]
###dmp <- GRanges("chr11",IRanges(tmp$cpg_pos,tmp$cpg_pos))
###dmpTrack <- AnnotationTrack(dmp,name="DMP")
###displayPars(dmpTrack) <- list(col=NULL,fill="darkblue")
###
###tmp <- unique(snpLoc$pos)
###snp <- GRanges("chr11", IRanges(tmp,tmp))
###snpTrack <- AnnotationTrack(snp,name="SNP")
###displayPars(snpTrack) <- list(col=NULL,fill="gray20")
###
###itrack <- IdeogramTrack(genome="hg19")
###gtrack <- GenomeAxisTrack(genome="hg19",cex=0.5)
###pdf(sprintf("%s/cis_eQTL_IGF2_Gviz.pdf",dirname(resFile)),width=6,height=3)
###plotTracks(list(itrack,intTrack,dmpTrack,snpTrack,refGenes,gtrack), 
###	chromosome="chr11", from=2100000,to=2200000,
###	sizes=list(0.1,0.3,0.1,0.1,0.3,0.1))
###plotTracks(list(itrack,intTrack,dmpTrack,snpTrack,refGenes,gtrack), 
###	chromosome="chr11", from=2150000,to=2200000,
###	sizes=list(0.1,0.3,0.1,0.1,0.3,0.1))
###dev.off()



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


