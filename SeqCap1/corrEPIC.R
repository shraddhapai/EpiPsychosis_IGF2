#' correlate methylation with EPIC
rm(list=ls())
require(minfi)

# both need to be in hg19 or hg38
minBaseCount <- 10
epicFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata"
seqCapDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCapEpi/methylation"
dmpFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/comb-p_171204/agesexPMI_PC12_171129_genes.txt"
minCvg <- 10 # min cvg for base to eligible for comparison

outDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCapEpi/DMR"

LIFTOVER <- "/home/shraddhapai/software/kent_utilities/liftOver"
LO_CHAIN <- "/home/shraddhapai/software/kent_utilities/hg19ToHg38.over.chain.gz"

dt <- format(Sys.Date(),"%y%m%d")
logFile <- sprintf("%s/CorrIGF_EPIC_%s.log",outDir,dt)
sink(logFile,split=TRUE)
tryCatch({

# must be in hg19
dmp <- read.delim(dmpFile,sep="\t",h=T,as.is=T)
dmp <- dmp[which(dmp$z_sidak_p< 0.05 & dmp[,1]=="chr11" & dmp$start < 2154500),]
#targetGR <- GRanges("chr11",IRanges(2016393,2160000))
targetGR <- GRanges(dmp[,1],IRanges(dmp[,2],dmp[,3]),
	name=sprintf("target%i",1:nrow(dmp)))

# --------------------------------------------------
# prepare EPIC data
# read methylation
load(epicFile)
locs <- getLocations(MSet.genome)
idx <- as.data.frame(findOverlaps(locs, targetGR))
cat(sprintf("Target region: %i probes\n", nrow(idx)))

MSet.genome <- MSet.genome[idx$queryHits,]
locs <- locs[idx$queryHits]
pd <- pData(MSet.genome)
idx <- which(!duplicated(pd$ID))
MSet.genome <- MSet.genome[,idx]
pd <- pData(MSet.genome)
betas <- getBeta(MSet.genome)

locs2 <- getLocations(MSet.genome)
idx <- as.data.frame(findOverlaps(locs2,targetGR))

# average by target region
betas <- aggregate(betas[idx$queryHits,],by=list(target=targetGR$name[idx$subjectHits]),FUN=mean)

require(reshape2)
betas <- melt(betas)
colnames(betas)[2] <- "Var2"
pd$Var2 <- rownames(pd)
pd <- pd[,c("ID","Sample_ID","DX","Var2")]
x <- merge(x=pd,y=betas, by="Var2")

locs <- as.data.frame(targetGR)
cat("Convert locs and targetGR to hg38\n")
locs$start <- locs$start-1 # open-1 position for ucsc
write.table(locs[,c("seqnames","start","end","name")],
	file="epic.hg19.txt",sep="\t",col=F,row=F,quote=F)
cmd <- sprintf("%s epic.hg19.txt %s epic.hg38.txt epic.unmapped.txt",
	LIFTOVER,LO_CHAIN)
system(cmd)
locs_hg38 <- read.delim("epic.hg38.txt",sep="\t",h=F,as.is=T)
colnames(locs_hg38) <- c("seqnames","start","end","Var1")
cat(sprintf("Loci: %i in hg19 -> %i in hg38", 
	nrow(locs),nrow(locs_hg38)))
targetGR <- GRanges(locs_hg38[,1],IRanges(locs_hg38[,2],locs_hg38[,3]),
	name=locs_hg38[,4])
rm(locs)

# now do same for target_GR
###tmp <- as.data.frame(targetGR)
###browser()
###write.table(tmp[,c(1:3)],file="target.hg19.txt",
###	sep="\t",col=F,row=F,quote=F)
###cmd <- sprintf("%s target.hg19.txt %s target.hg38.txt target.unmapped.txt",
###	LIFTOVER,LO_CHAIN)
###system(cmd)
###target.hg19 <- tmp
###tmp <- read.delim("target.hg38.txt",sep="\t",h=F,as.is=T)
###targetGR <- GRanges(tmp[,1],IRanges(tmp[,2],tmp[,3]))

x$Var1 <- as.character(x$target)
y <- merge(x=x,y=locs_hg38,by="Var1")
y <- as.data.frame(y) # convert from DataFrame object

require(ggplot2)
p <- ggplot(y,aes(x=start,y=value))+geom_point(aes(colour=DX),cex=0.1) + geom_smooth(aes(colour=DX))+geom_hline(yintercept=0.5,lty=2)
epic_vals <- y

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
for (fName in fList) {
	sampName <- sub(".CaseControl.records.Rdata","",fName)
	print(sampName)
	out <- poolStrands(sprintf("%s/%s",seqCapDir,fName),getTargetGR=TRUE)
	# initial filter to bases in targets
	idx <- as.data.frame(findOverlaps(targetGR,out$target_GR))
	rec2 <- do.call("rbind",out$rec[idx$subjectHits])
	rec2 <- rec2[!duplicated(rec2),]
	rec2 <- subset(rec2[which(rec2$CT_count>=minBaseCount),])
	rec_GR <- GRanges(rec2[,1],IRanges(rec2[,2],rec2[,2]))
	rec_GR$C_count <- rec2$C_count
	rec_GR$CT_count <- rec2$CT_count
	rec2$pctM <- rec2$C_count/rec2$CT_count
	
	# now average by targets
	ol <- as.data.frame(findOverlaps(rec_GR,targetGR))
	mu <- aggregate(rec2[ol$queryHits,c("C_count","CT_count","pctM")],
		by=list(target=targetGR$name[ol$subjectHits]),
		FUN=mean)
	tmp <- melt(mu)
	rsamp <- sub("-[123]re","",sampName)
	tmp$Sample_ID <- rsamp
	seqcap_vals[[fName]] <- tmp
}

# data.frame should have region, sample, SeqCap, EPIC
# output graph should show x-axis EPIC, y-axis SeqCap
# each dot should be one sample x one base

epic_vals <- epic_vals[,c("target","Sample_ID","value","DX")]
colnames(epic_vals)[3] <- "avgM_EPIC"
seqcap_vals <- do.call("rbind",seqcap_vals)
seqcap_vals$Sample_ID <- sub("re","",seqcap_vals$Sample_ID)
colnames(seqcap_vals)[3] <- "avgPctM_SeqCapEPI"
seqcap_vals <- subset(seqcap_vals, variable=="pctM")
agg <- aggregate(seqcap_vals$avgPctM_SeqCapEPI,
	by=list(Sample_ID=seqcap_vals$Sample_ID, target=seqcap_vals$target),
	FUN=mean)
colnames(agg)[3] <- "avgPctM_SeqCapEPI"

#browser()

comb <- merge(x=epic_vals,y=seqcap_vals,by=c("target","Sample_ID"))
comb <- subset(comb, variable =="pctM")
comb$Sample_ID <- factor(comb$Sample_ID)
# mixed effects model
#fit0 <- lmer(avgM_EPIC~1+(1|Sample_ID),data=comb,REML=TRUE)
#fit1 <- lmer(avgM_EPIC~1+avgPctM_SeqCapEPI+(1|Sample_ID),data=comb,REML=TRUE)
#x <- anova(fit1,fit0)
#cat(sprintf("LME: ANOVA p < %1.2e\n", x[["Pr(>Chisq)"]][2]))


comb <- merge(x=epic_vals,y=agg,by=c("target","Sample_ID"))
comb$DX <- factor(comb$DX,levels=c("case","control"))
fit <- lm(avgM_EPIC~avgPctM_SeqCapEPI,data=comb)
print(summary(fit))

comb$avgPctM_SeqCapEPI <- comb$avgPctM_SeqCapEPI
p <- ggplot(comb, aes(x=avgPctM_SeqCapEPI,y=avgM_EPIC))
p <- p + geom_point(aes(colour=DX),cex=4) + ylim(c(0,1)) + xlim(c(0,1))
p <- p + geom_hline(yintercept=0,lty=2,lwd=2) 
p <- p + geom_vline(xintercept=0,lty=2,lwd=2)
p <- p + geom_abline(slope=1,intercept=0,lwd=1.5)
p <- p + geom_abline(slope=coef(fit)[2], intercept=coef(fit)[1],
	lty=2,col='grey50',lwd=1.5)
p <- p + ggtitle(sprintf("IGF2: %s:%i-%i\nEPIC = %1.2f * SeqCap + %1.2f",
		seqnames(targetGR), min(start(targetGR)), max(end(targetGR)),
		coef(fit)[2],coef(fit)[1]))
p <- p + theme(text=element_text(size=20),legend.text=element_text(size=20))
pdf(sprintf("%s/corrIGF2_DMR_%s.pdf",outDir,dt));print(p); dev.off()

},error=function(ex){ print(ex) }, finally={sink(NULL)})
