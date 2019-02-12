# volcano plot for methylation analysis
rm(list=ls())

require(dataExplore) # qqplot

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
diffExFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlRNAseq/SP_Lee_compare/STAR_SP_output/STAR_edgeR_GroupAgeSexPMINeurons_180111.txt"

outDir <- dirname(diffExFile)
dmp <- read.delim(diffExFile,sep="\t",h=T,as.is=T)
colnames(dmp)[which(colnames(dmp)=="V6")] <- "gene.name"

dt <- format(Sys.Date(),"%y%m%d")

postscript(sprintf("%s/dmp_qqplot_%s.eps", outDir,dt))
ggd.qqplot(dmp$PValue,main="nominal p, qqplot")
dev.off()
png(sprintf("%s/dmp_qqplot_%s.png", outDir,dt))
ggd.qqplot(dmp$PValue,main="nominal p, qqplot",cex.axis=2)
dev.off()

postscript(sprintf("%s/dmp_pvaluehistogram_%s.eps", outDir,dt))
hist(dmp$PValue,n=100,main="nominal p histogram")
dev.off()

dmp$P <- -log10(dmp$PValue)
idx <- which(dmp$P >4)
dmp2 <- dmp[idx,]
dmp2$names <- dmp2$gene.name#rownames(dmp2)

dmp3 <- dmp[which(dmp$gene.name %in% "IGF2"),,drop=FALSE]

require(ggplot2)
# volcano plot
cat("* Creating volcano plot\n")
p <- ggplot(dmp,aes(x=logFC,y=P))+geom_point(size=0.3,colour="grey50")
p <- p + geom_point(data=dmp2,aes(x=logFC,y=P),colour="red",size=1.5)
p <- p + theme(text=element_text(size=24),axis.text=element_text(size=24))
p <- p + ggtitle(sprintf("Psychosis gene diffEx (N=%s genes)",
		prettyNum(nrow(dmp),big.mark=",")))
p <- p +  ggrepel::geom_text_repel(
    data = dmp2,
    aes(label = names),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
# add IGF2
p <- p + geom_point(data=dmp3,aes(x=logFC,y=P),colour="blue",size=2)
p <- p +  ggrepel::geom_text_repel(
    data = dmp3,
    aes(label = gene.name),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )



cat("* Writing to file\n")
postscript(sprintf("%s/dmp_volcanoPlot_%s.eps",outDir,dt),width=6,height=6)
print(p)
dev.off()
png(sprintf("%s/dmp_volcanoPlot_%s.png",outDir,dt))
print(p)
dev.off()

# histogram of fold-change
cat("Writing fc histogram\n")
dmp3 <- subset(dmp,PValue<0.05)
dmp3$FC <- 2^dmp3$logFC
p <- ggplot(dmp3,aes(x=FC)) + geom_histogram(colour="white")
p <- p + theme(text=element_text(size=16),axis.text=element_text(size=14))
p <- p + ggtitle(sprintf("Fold change, probes with P < 0.05 (N=%i)", 
	nrow(dmp3)))
p <- p + geom_vline(xintercept=1,lty=2)
postscript(sprintf("%s/dmp_FChistogram_%s.eps",outDir,dt),width=6,height=6)
print(p); dev.off()
png(sprintf("%s/dmp_FChistogram_%s.png",outDir,dt))
print(p); dev.off()




