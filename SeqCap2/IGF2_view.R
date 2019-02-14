#' look at IGF2 and TH regions.
require(GenomicRanges)
source("poolStrands.R")
require(reshape2)
require(limma)
require(ggplot2)

inDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/methylation/locuswise"
outDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/DMR"
phenoFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/SeqCap2_pheno_v3.txt"
phenoTechFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/SeqCap2_samplePheno.txt"
minCvg <- 10

excSamples <- TRUE
dt <- format(Sys.Date(),"%y%m%d")
logFile <- sprintf("%s/IGF2_DMR_view_exc%s_%s.log",outDir,excSamples,dt)
sink(logFile,split=TRUE)

tryCatch({
fList <- dir(inDir, pattern="CaseControl.records.Rdata")
methList <- list()
for (fName in fList) {
	sampName <- sub(".CaseControl.records.Rdata","",fName)
	print(sampName)
	out <- poolStrands(sprintf("%s/%s",inDir,fName),getTargetGR=TRUE)
	rec <- out$rec
	# target-level methylation
	pctM <- lapply(rec,function(x) {
		if (!is.null(dim(x))) {
			idx <- which(x$CT_count >= minCvg)	
			if (any(idx)) {
				y <- sum(x$C_count[idx])/sum(x$CT_count[idx])
				return(y)
			}
			else return(NA)
		} else return(NA)
	})
	pctM <- unlist(pctM)
	names(pctM) <- names(rec)
	methList[[sampName]] <- pctM	
}

	# compile into a table - rows are regions, columns are samples
	# cells are overall pctM
	uq_rg <- names(methList[[1]])
	rgn_meth <- matrix(NA,nrow=length(uq_rg), ncol=length(methList))
	rownames(rgn_meth) <- uq_rg
	cat("Compiling table\n")
	for (k in 1:length(methList)) {
		cur <- methList[[k]]
		midx <- match(names(cur),rownames(rgn_meth))
		if (all.equal(names(cur),rownames(rgn_meth)[midx])!=TRUE) {
			cat("compile table, don't match\n")
		}
		rgn_meth[midx,k] <- cur
	}
	colnames(rgn_meth) <- names(methList)

# run hierarchical clustering of case-control loci.
x2 <- na.omit(rgn_meth)
x <- as.dist(1-cor(x2))
pdf(sprintf("%s/hclust_exc%s_%s.pdf",outDir,excSamples,dt)) 
plot(hclust(x,method="average"),
	main=sprintf("SeqCap2: %i DMR candidate loci",nrow(x2)))
dev.off()

# -------------------
# colourful hclust plot
pheno <- read.delim(phenoTechFile,sep="\t",h=T,as.is=T)

colnames(rgn_meth)[which(colnames(rgn_meth)=="95repos")] <- "90-2repos"
sampName <- sub("Redo","",colnames(rgn_meth))
sampName <- sub("pos","",sampName)
sampName <- sub("neg","",sampName)
sampName <- sub("-[123]","",sampName)
sampName <- as.integer(sub("re","",sampName))
pheno$New.Code <- pheno$Patient.ID 
midx <- match(sampName,pheno$New.Code)
if (all.equal(pheno$New.Code[midx],sampName)!=TRUE) {
	cat("pheno sampName order doesn't match"); browser()
}
pdf(sprintf("%s/plotDendro_before_%s.pdf",outDir,dt))
dev.off()

if (excSamples) {
	idx <- which(colnames(rgn_meth) %in% c("75pos","92pos","9Redopos"))
	rgn_meth <- rgn_meth[,-idx]
	sampName <- sampName[-idx]
}
midx <- match(sampName, pheno$New.Code)
if (all.equal(pheno$New.Code[midx],sampName)!=TRUE) {
	cat("pheno sampName order doesn't match"); browser()
}

x2 <- na.omit(rgn_meth)
x <- as.dist(1-cor(x2))
pdf(sprintf("%s/hclust_exc%s_after_removing_%s.pdf",outDir,excSamples,dt))
plot(hclust(x,method="average"),
    main=sprintf("SeqCap2: %i DMR candidate loci",nrow(x2)))
dev.off()


# ----------------------------------------------------
# Compare IGF2 in neurons vs glia
#90 and 95 are technical replicates

df <- melt(rgn_meth)
df$value <- df$value*100

pheno <- read.delim(phenoFile)
sampName <- sub("pos","",df$Var2)
sampName <- sub("neg", "",sampName)
cellType <- rep("",nrow(df))
cellType[grep("neg",df$Var2)] <- "glia"
cellType[grep("pos",df$Var2)] <- "neuron"
df$cellType <- cellType
df$ID <- sampName
df$Var2 <- sub("-[123]re","",df$Var2)
df$Var2 <- sub("re","",df$Var2)
df$ID <- sub("-[123]re","",df$ID)
df$ID <- sub("re","",df$ID)

agg <- aggregate(df$value,by=list(ID=df$ID,cellType=df$cellType,Var1=df$Var1),
	FUN=mean,na.rm=T)
colnames(agg)[4] <- "value"
df <- agg

colnames(pheno)[1] <- "ID"
y <- merge(x=pheno,y=df,by="ID")
y <- y[!duplicated(y),]
ph <- y

ph$cellType <- factor(ph$cellType)
ph$SEX <- factor(ph$SEX)
ph$New.Code <- factor(ph$New.Code)
dx <- rep(NA,nrow(ph))
dx[which(ph$DIST.DX %in% "Control")] <- "control"
dx[which(ph$DIST.DX %in% c("Schizophrenia","Bipolar","Bipolar Disorder"))] <- "case"
ph$DX <- factor(dx)

ph <- subset(ph, Var1 %in% "CaseControl_12_3")
neuron <- subset(ph,cellType=="neuron")
colnames(neuron)[which(colnames(neuron) == "value")] <- "neuron"
glia <- subset(ph,cellType=="glia")[,c("New.Code","value")]
colnames(glia)[which(colnames(glia) == "value")] <- "glia"
both <- merge(x=neuron,y=glia,by=c("New.Code"),all.x=TRUE,all.y=TRUE)

missid <- as.character(unique(both$New.Code[which(is.na(both$SEX))]))
for (k in missid) {
	idx <- which(both$New.Code==k)
	orig <- ph[which(pheno$New.Code == k)[1],]
	both$SEX[idx] <- orig$SEX
	both$RACE[idx] <- orig$RACE
	both$AGE[idx] <- orig$AGE
	both$PMI[idx] <- orig$PMI
	both$ID[idx] <- orig$ID
	both$batch[idx] <- orig$batch
}

both <- both[!duplicated(both),]
both$cellDiff <- both$neuron-both$glia
both$DX <- factor(both$DX)

updatedSampKey <- "/home/shraddhapai/Epigenetics/NARSAD/input_files/NARSAD_sampleKey_171129.txt"

#### update PMI for last three samples
p2 <- read.delim(updatedSampKey,sep="\t",h=T,as.is=T)
both$PMI <- as.numeric(as.character(both$PMI))
idx <- which(is.na(both$PMI)) # find missing pmi
cat("Filling in missing PMIs\n")
for (k in idx) {
	idx2 <- which(as.character(p2$Internal.ID) == both$New.Code[k])
    cat(sprintf("\tSample %s\n", both$New.Code[k]))
	both$PMI[k] <- p2$PMI[idx2]
}

require(ggplot2)
require(GGally)
both2 <- subset(both, Var1 %in% "CaseControl_12_3")
p <- ggplot(both2,aes(x=DX,y=cellDiff))
p <- p + geom_boxplot(aes(fill=SEX),lwd=1.5,cex.axis=1.5)
p <- p + geom_hline(yintercept=0) + ylab("Neuron - Glia")
p <- p + ggtitle("IGF2 - celltype difference")
p <- p + theme_bw() + theme(axis.text=element_text(size=16))
pdf(sprintf("%s/IGF2_CellTypeDiff_%s_%s.pdf",outDir,excSamples,dt)); 
ggparcoord(both2,columns=c(which(colnames(both)%in% "neuron"), 
	which(colnames(both)%in%"glia")),
	groupColumn="New.Code",scale="globalminmax")
print(p); 
dev.off()

cat("Testing Neuron-Glia for IGF2\n")
print(table(both2$DIST.DX))
y <- t.test(both2$cellDiff,alternative="greater")
print(y)


# now look at IGF2 case/control in neurons, vs that in glia
ph2 <- subset(ph, Var1 %in% "CaseControl_12_3")
p <- ggplot(ph2,aes(x=cellType,y=value))
p <- p + geom_boxplot(aes(fill=DX),lwd=1.5)
p <- p + ggtitle("case/control IGF2: split by celltype") 
p <- p + ylim(0,55)
p <- p + theme_bw() + theme(axis.text=element_text(size=16))
pdf(sprintf("%s/IGF2_CaseControl_NeuronGlia_%s_%s.pdf",outDir,excSamples,dt)); 
print(p); dev.off()

ph2$DIST.DX[which(ph2$DIST.DX %in% "Bipolar")] <- "Bipolar Disorder"

#### update PMI for last three samples
p2 <- read.delim(updatedSampKey,sep="\t",h=T,as.is=T)
ph2$PMI <- as.numeric(as.character(ph2$PMI))
idx <- which(is.na(ph2$PMI)) # find missing pmi
cat("Filling in missing PMIs\n")
for (k in idx) {
	idx2 <- which(as.character(p2$Internal.ID) == ph2$New.Code[k])
    cat(sprintf("\tSample %s\n", ph2$New.Code[k]))
	ph2$PMI[k] <- p2$PMI[idx2]
}

cat("--------\n")
# test IGF2 enhancer in neuron
neudat <- subset(ph2, cellType %in% "neuron")
cat(sprintf("\tNeuron: %i data points\n", nrow(neudat)))
print(table(neudat$DX))
print(table(neudat$DIST.DX))
fit0 <- lm(value~1+AGE+SEX+PMI+Batch,data=neudat)
fit1 <- lm(value~1+DX+AGE+SEX+PMI+Batch,data=neudat)
res <- anova(fit1,fit0)
neu_anova <- res[["Pr(>F)"]][2]
browser()
cat(sprintf("ANOVA: p < %1.2e\n", res[["Pr(>F)"]][2]))
neu_mudiff <- mean(neudat$value[which(neudat$DX %in% "case")])-mean(neudat$value[which(neudat$DX %in% "control")])
cat(sprintf("Mean case - Mean control = %1.2f %%\n", neu_mudiff))

cat("--------\n")
gldat <- subset(ph2,cellType %in% "glia")
cat(sprintf("\tGlia: %i data points\n", nrow(gldat)))
print(table(gldat$DX))
print(table(gldat$DIST.DX))
fit0 <- lm(value~1+AGE+SEX+PMI,data=gldat)
fit1 <- lm(value~1+DX+AGE+SEX+PMI,data=gldat)
res <- anova(fit1,fit0)
glia_anova <- res[["Pr(>F)"]][2]
cat(sprintf("ANOVA: p < %1.2e\n", res[["Pr(>F)"]][2]))
gl_mudiff <- mean(gldat$value[which(gldat$DX %in% "case")])-mean(gldat$value[which(gldat$DX %in% "control")])
cat(sprintf("Mean case - Mean control = %1.2f %%\n", gl_mudiff))

}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})

