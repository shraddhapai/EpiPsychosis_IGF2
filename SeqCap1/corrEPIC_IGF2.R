#' correlate methylation with EPIC
rm(list=ls())
require(minfi)

# both need to be in hg19 or hg38
minBaseCount <- 10
epicFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata"
seqCapDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCapEpi/methylation"
outDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCapEpi/DMR"
LIFTOVER <- "/home/shraddhapai/software/kent_utilities/liftOver"
LO_CHAIN <- "/home/shraddhapai/software/kent_utilities/hg19ToHg38.over.chain.gz"

dt <- format(Sys.Date(),"%y%m%d")
logFile <- sprintf("%s/CorrIGFextended_EPIC_%s.log",outDir,dt)
sink(logFile,split=TRUE)
tryCatch({
# --------------------------------------------------
# now getting SeqCap data

#### get sample-level methylation data for each base
#### get it for exact base? merge nearby bases?
#### convert to table
#### now correlate with EPIC
source("getM_GRanges.R")
source("poolStrands.R")

# 31-3re moved out of this folder because it doesn't cluster with its
# technical replicates
fList <- dir(seqCapDir, pattern="CaseControl.records.Rdata")

seqcap_vals <- list()
seqcap_range <- NULL
for (fName in fList) {
	sampName <- sub(".CaseControl.records.Rdata","",fName)
	print(sampName)
	out <- poolStrands(sprintf("%s/%s",seqCapDir,fName),getTargetGR=TRUE)
	# initial filter to bases in targets
	idx <- which(out$target_GR$name=="CaseControl_12_3")
	if (is.null(seqcap_range)) {
		seqcap_range <- out$target_GR[idx]
	}
	rec2 <- out$rec[[idx]]
	rec2 <- rec2[!duplicated(rec2),]
	rec2 <- subset(rec2[which(rec2$CT_count>=minBaseCount),])
	rec2$pctM <- rec2$C_count/rec2$CT_count
	
	rsamp <- sub("-[123]re","",sampName)
	rsamp <- sub("re","",rsamp)
	rec2$Sample_ID <- rsamp
	rec2$start <- rec2$pos
	rec2$pos <- sprintf("%s:%i-%i",rec2$chr,rec2$pos-1,rec2$pos-1)
	seqcap_vals[[fName]] <- rec2[,c("pos","Sample_ID","pctM","start")]
}
seqcap_vals <- do.call("rbind",seqcap_vals)
# aggregate by tech reps
agg <- aggregate(seqcap_vals$pctM, 
	by=list(Sample_ID=seqcap_vals$Sample_ID,pos=seqcap_vals$pos),
	FUN=mean)
colnames(agg)[3] <- "pctM"
seqcap_vals <- agg

rownames(seqcap_vals) <- NULL
rm(rec2)

# --------------------------------------------------
# prepare EPIC data
# read methylation
load(epicFile)
locs <- getLocations(MSet.genome)
hg19_zone <- GRanges("chr11",IRanges(2147342,2165341))

idx <- as.data.frame(findOverlaps(locs, hg19_zone))
cat(sprintf("Target region: %i probes\n", nrow(idx)))

MSet.genome <- MSet.genome[idx$queryHits,]
locs <- as.data.frame(getLocations(MSet.genome))
# lift over to hg38
locs$start <- locs$start-1 # open-1 position for ucsc
locs$name <- rownames(locs)
write.table(locs[,c("seqnames","start","end","name")],
	file="epic.hg19.txt",sep="\t",col=F,row=F,quote=F)
cmd <- sprintf("%s epic.hg19.txt %s epic.hg38.txt epic.unmapped.txt",
	LIFTOVER,LO_CHAIN)
system(cmd)
locs_hg38 <- read.delim("epic.hg38.txt",sep="\t",h=F,as.is=T)
colnames(locs_hg38) <- c("seqnames","start","end","Var1")
cat(sprintf("Loci: %i in hg19 -> %i in hg38", nrow(locs),nrow(locs_hg38)))

if (all.equal(locs$name,locs_hg38[,4])!=TRUE) {
	cat("ids don't match\n")
	browser()
}

locs <- locs_hg38
locs_pos <- sprintf("%s:%i-%i",locs$seqnames,locs$start,locs$start)
locs$pos <- locs_pos
betas <- getBeta(MSet.genome)
pd <- pData(MSet.genome)
colnames(betas) <- pd$Sample_ID
if (all.equal(rownames(betas),locs$Var1)!=TRUE) {
	cat("beta locs don't match"); browser()
}

betas <- cbind(locs,betas)
betas <- betas[,-(1:3)]
betas2 <- melt(betas)
colnames(betas2)[3] <- "Sample_ID"

# plot view to confirm how seqcap and epic coordinates align
# e.g. is one of them off by 1 consistently, rel to other?
###uq <- unique(seqcap_vals$start)
###xmin <- 2132900; xmax <- 2133000
###pdf("test.pdf",width=11,height=4)
###plot(locs$start,rep(1,nrow(locs)),pch=16,cex=0.5,ylim=c(0,3),xlim=c(xmin,xmax))
###points(uq,rep(2,length(uq)),pch=16,cex=0.5,col='red')
###segments(x0=xmin:xmax,y0=0,y1=3,lty=3,col='grey50')
###dev.off()

y <- merge(x=betas2,y=seqcap_vals,by=c("pos","Sample_ID"))
cor_p <- cor.test(y$value,y$pctM,method="p")
cor_sp <- cor(y$value,y$pctM,method="sp")
cat(sprintf("Correlation is %1.2f Pearson, %1.2f Spearman\n",cor_p$estimate,cor_sp))

p <- ggplot(y, aes(x=pctM,y=value)) + geom_point(aes(colour=Sample_ID)) 
p <- p+ ylab("M value from epic") + xlab("seqcap pctM") 
p <- p + xlim(c(0,1)) + ylim(c(0,1))
p <- p + geom_abline(slope=1,intercept=0) 
p <- p + ggtitle(sprintf("IGF2,seqcap/epic (minCvg>=%i)\nCorr: p:%1.2f (p < %1.2e)",
	minBaseCount,cor_p$estimate,cor_p$p.value))
outFile <- sprintf("%s/corrSeqCapEPI_IGF2_minCvg%i_%s.pdf",
	outDir,minBaseCount,dt)
pdf(outFile) print(p); dev.off()


},error=function(ex){print(ex)},finally={ sink(NULL)})
