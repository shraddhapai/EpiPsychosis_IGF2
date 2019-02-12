# volcano plot for methylation analysis
rm(list=ls())
require(minfi)
#dmpFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/agesexslide_PC12/161221/dmp..dmp_DXcase.topTable.161221.txt"
dmpFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/agesexPMI_PC12_171129/dmp..dmp_DXcase.topTable.171129.txt"
combPFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/out_agesexPMI_PC12_171129.regions-p.bed"
resFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata"

# -----------------------------------------
# Utility functions
# Define the function
ggd.qqplot = function(pvector, main=NULL, ...) {
		mar <- par("mar")
		par(mar=c(5,6,2,1))
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o)),
				cex.axis=3,cex.lab=2)
    lines(e,e,col="red")
		par(mar=mar)
}

# -----------------------------------------


outDir <- dirname(dmpFile)


combP <- read.delim(combPFile,sep="\t",h=T,as.is=T)
combP <- subset(combP, z_sidak_p < 0.05)
require(GenomicRanges)
combGR <- GRanges(combP[,1],IRanges(combP[,2],combP[,3]))

dmp <- read.delim(dmpFile,sep="\t",h=T,as.is=T)

load(resFile)
mLocs <- getLocations(MSet.genome)
redP <- subsetByOverlaps(mLocs,combGR)
redP <- names(redP)

dmp$P <- -log10(dmp$P.Value)

# outdated
# idx <- which(dmp$adj.P.Val < 0.05 & rownames(dmp) %in% probeIDs)
idx <- which(rownames(dmp) %in% redP)
dmp2 <- dmp[idx,]
dmp2$names <- rownames(dmp2)


probeIDs <- c("cg07096953","cg26401390","cg02613624","cg22956483")
dmp3 <- dmp[which(rownames(dmp) %in% probeIDs),]

require(ggplot2)

# volcano plot
cat("* Creating volcano plot\n")
p <- ggplot(dmp,aes(x=logFC,y=P))+geom_point(size=0.3,colour="grey50")
# probes in comb-p sig regions are red
p <- p + geom_point(data=dmp2,aes(x=logFC,y=P),colour="red",size=1.5)
p <- p + geom_point(data=dmp3,aes(x=logFC,y=P),colour="blue",size=2.5)
p <- p + theme(text=element_text(size=24),axis.text=element_text(size=24))
p <- p + ggtitle(sprintf("Psychosis methylation changes (N=%s probes)",
		prettyNum(nrow(dmp),big.mark=",")))
#p <- p +  ggrepel::geom_text_repel(
#    data = dmp2,
#    aes(label = names),
#    size = 5,
#    box.padding = unit(0.35, "lines"),
#    point.padding = unit(0.3, "lines")
#  )
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

cat("* Writing to file\n")
dt <- format(Sys.Date(),"%y%m%d")
postscript(sprintf("%s/dmp_volcanoPlot_%s.eps",outDir,dt),width=6,height=6)
print(p)
dev.off()
png(sprintf("%s/dmp_volcanoPlot_%s.png",outDir,dt))
print(p)
dev.off()

# qqplot
postscript(sprintf("%s/dmp_qqplot_%s.eps", outDir,dt))
ggd.qqplot(dmp$P.Value,main="nominal p, qqplot")
dev.off()
png(sprintf("%s/dmp_qqplot_%s.png", outDir,dt))
ggd.qqplot(dmp$P.Value,main="nominal p, qqplot",cex.axis=2)
dev.off()


# histogram of fold-change
cat("Writing fc histogram\n")
dmp3 <- subset(dmp,P.Value<0.05)
dmp3$FC <- 2^dmp3$logFC
p <- ggplot(dmp3,aes(x=FC)) + geom_histogram(colour="white")
p <- p + theme(text=element_text(size=8),axis.text=element_text(size=14))
p <- p + ggtitle(sprintf("Fold change, probes with P < 0.05 (N=%i)", 
	nrow(dmp3)))
p <- p + geom_vline(xintercept=1,lty=2)
postscript(sprintf("%s/dmp_FChistogram_%s.eps",outDir,dt),width=6,height=6)
print(p); dev.off()
png(sprintf("%s/dmp_FChistogram_%s.png",outDir,dt))
print(p); dev.off()




