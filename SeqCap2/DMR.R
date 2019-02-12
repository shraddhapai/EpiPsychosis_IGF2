#' look at IGF2 and TH regions.
rm(list=ls())
require(GenomicRanges)
require(reshape2)
require(limma)
require(ggplot2)

inDir <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/methylation/locuswise"
phenoFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/SeqCap2_pheno_v3.txt"
phenoTechFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/SeqCap2_samplePheno.txt"
updatedSampKey <- "/home/shraddhapai/Epigenetics/NARSAD/input_files/NARSAD_sampleKey_171129.txt"

runDMR <- function(gpName,cellType_in) {
source("poolStrands.R")
minCvg <- 10
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/DMR/%s_%s_%s",
	gpName,cellType_in,dt)
if (!file.exists(outDir)) dir.create(outDir)

excSamples <- TRUE
logFile <- sprintf("%s/DMR_%s_view_exc%s_%s.log",outDir,gpName,excSamples,dt)
sink(logFile,split=TRUE)

tryCatch({
fList <- dir(inDir, pattern=sprintf("%s.records.Rdata",gpName))
cat(sprintf("Group: %s -- got %i files\n", gpName, length(fList)))
methList <- list()
wroteTargets <- FALSE
tgtdf <- NULL
for (fName in fList) {
	sampName <- sub(sprintf(".%s.records.Rdata",gpName),"",fName)
	print(sampName)
	out <- poolStrands(sprintf("%s/%s",inDir,fName),getTargetGR=TRUE)
	if (!wroteTargets) {
		tgtdf <- as.data.frame(out$target_GR)
		write.table(tgtdf,file=sprintf("%s/target_gr.txt",outDir),sep="\t",
			col=T,row=F,quote=F)
		wroteTargets <- TRUE
	}
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

write.table(rgn_meth,file=sprintf("%s/methylation.txt",outDir),sep="\t",col=T,row=T,quote=F)


# run hierarchical clustering of case-control loci.
x2 <- na.omit(rgn_meth)
###pdf(sprintf("%s/hclust_exc%s_%s.pdf",outDir,excSamples,dt)) 
###plot(hclust(x,method="average"),
###	main=sprintf("SeqCap2: %s: %i DMR candidate loci",gpName,nrow(x2)))
###dev.off()

# -------------------
# colourful hclust plot
pheno <- read.delim(phenoTechFile,sep="\t",h=T,as.is=T)
#browser()

colnames(rgn_meth)[which(colnames(rgn_meth)=="95repos")] <- "90-2repos"
sampName <- sub("Redo","",colnames(rgn_meth))
sampName <- sub("pos","",sampName)
sampName <- sub("neg","",sampName)
sampName <- sub("-[123]","",sampName)
sampName <- as.integer(sub("re","",sampName))
midx <- match(sampName,pheno$Patient.ID)
if (all.equal(pheno$Patient.ID[midx],sampName)!=TRUE) {
	cat("pheno sampName order doesn't match"); browser()
}
pheno <- pheno[midx,]
pdf(sprintf("%s/plotDendro_before_%s.pdf",outDir,dt))
x <- as.dist(1-cor(x2))
plot(hclust(x,method="average"),
    main=sprintf("SeqCap2: %i DMR candidate loci",nrow(x2)))
dev.off()
#plotDendro_clr(na.omit(rgn_meth),pheno[midx],

if (excSamples) {
	idx <- which(colnames(rgn_meth) %in% c("75pos","92pos","9Redopos"))
	if (any(idx)) { ### check put in so debugging with fewer samples
					### doesn't break code
		rgn_meth <- rgn_meth[,-idx]
		sampName <- sampName[-idx]
	}
}
midx <- match(sampName,pheno$Patient.ID)
if (all.equal(pheno$Patient.ID[midx],sampName)!=TRUE) {
	cat("pheno sampName order doesn't match"); browser()
}
pheno <- pheno[midx,]


x2 <- na.omit(rgn_meth)
x <- as.dist(1-cor(x2))
pdf(sprintf("%s/hclust_exc%s_after_removing_%s.pdf",outDir,excSamples,dt))
plot(hclust(x,method="average"),
    main=sprintf("SeqCap2: %i DMR candidate loci",nrow(x2)))
dev.off()


# ----------------------------------------------------
#90 and 95 are technical replicates

df <- melt(rgn_meth)
df$value <- df$value*100

pheno <- read.delim(phenoFile)
sampName <- sub("pos","",df$Var2)
sampName <- sub("neg", "",sampName)
cellType <- rep("neuron",nrow(df)) # batch 1 were neurons 
cellType[grep("neg",df$Var2)] <- "glia"
cellType[grep("pos",df$Var2)] <- "neuron"
df$cellType <- cellType
df$ID <- sampName
df$Var2 <- sub("-[123]re","",df$Var2)
df$Var2 <- sub("re","",df$Var2)
df$ID <- sub("-[123]re","",df$ID)
df$ID <- sub("re","",df$ID)
print("markera")

cat(sprintf("Subsetting for cell type = %s\n",cellType))
df <- subset(df, cellType==cellType_in)
cat(sprintf("%i samples left\n", length(unique(df$ID))))

#if (length(duplicated(df[,c("ID","Var1","cellType")]))>0) {
#	cat("\tTaking average of technical replicates\n")
	agg <- aggregate(df$value,by=list(ID=df$ID,
		Var1=df$Var1),
		FUN=mean,na.rm=T)
	colnames(agg)[3] <- "value"
	df <- agg
#}

colnames(pheno)[1] <- "ID"
y <- merge(x=pheno,y=df,by="ID")
y <- y[!duplicated(y),]
ph <- y

ph$SEX <- factor(ph$SEX)
ph$New.Code <- factor(ph$New.Code)
dx <- rep(NA,nrow(ph))
dx[which(ph$DIST.DX %in% "Control")] <- "control"
dx[which(ph$DIST.DX %in% c("Schizophrenia","Bipolar","Bipolar Disorder"))] <- "case"
ph$DX <- factor(dx)
cat("markerb\n")

#### update PMI for last three samples
p2 <- read.delim(updatedSampKey,sep="\t",h=T,as.is=T)
ph$PMI <- as.numeric(as.character(ph$PMI))
idx <- which(is.na(ph$PMI)) # find missing pmi
cat("Filling in missing PMIs\n")
for (k in idx) {
	idx2 <- which(as.character(p2$Internal.ID) == ph$New.Code[k])
    cat(sprintf("\tSample %s\n", ph$New.Code[k]))
	ph$PMI[k] <- p2$PMI[idx2]
}

x <- dcast(data=ph,Var1~New.Code,value.var="value")
rownames(x) <- x[,1]; x <- x[,-1]
dat <- na.omit(x); rm(x)
var <- sort(dataExplore::getVariance(dat),decreasing=TRUE)
tokeep <- names(var)[1:(round(0.5*length(var)))]

cat(sprintf("Keeping only top 50%% most variable = %i loci\n",
	length(tokeep)))

ph <- subset(ph, Var1 %in% tokeep)
ph$DIST.DX[which(ph$DIST.DX %in% "Bipolar")] <- "Bipolar Disorder"
# locus-wise case/control in neurons, vs that in glia
out <- list()
numsamp <- length(unique(ph$New.Code))
for (tgt in unique(ph$Var1)) {
	print(tgt)
	neudat <- subset(ph, Var1 %in% tgt) # & cellType %in% "neuron")
	isnan <- which(is.na(neudat$value) | is.nan(neudat$value))
	if (any(isnan)) {
		neudat <- neudat[-isnan,]
	}
	if (nrow(neudat)>=0.75*numsamp) {
		cat(sprintf("\tNeuron: %i data points\n", nrow(neudat)))
		print(table(neudat$DX))
		fit0 <- lm(value~1+AGE+SEX+PMI+Batch,data=neudat)
		fit1 <- lm(value~1+DX+AGE+SEX+PMI+Batch,data=neudat)
		res <- anova(fit1,fit0)
		neu_anova <- res[["Pr(>F)"]][2]
		cat(sprintf("ANOVA: p < %1.2e\n", res[["Pr(>F)"]][2]))
		neu_mudiff <- mean(neudat$value[which(neudat$DX %in% "case")])-mean(neudat$value[which(neudat$DX %in% "control")])
		#cat(sprintf("Mean case - Mean control = %1.2f %%\n", neu_mudiff))
		
		# target, N tot, N case, N control, mean diff, pvalue
		out[[tgt]] <- c(tgt,nrow(neudat), 
			sum(neudat$DX %in% "case"),sum(neudat$DX %in% "control"), 
			neu_mudiff,neu_anova)	
	}
}
out2 <- as.data.frame(do.call("rbind",out))
out2[,6] <- as.numeric(as.character(out2[,6]))
pdfFile <- sprintf("%s/nominalp.pdf",outDir)
pdf(pdfFile)
hist(out2[,6],n=100);dev.off()

out2$qval <- p.adjust(out2[,6],method="BH")
colnames(out2) <- c("target","N","N_case","N_ctrl","mean(case-ctrl)",
	"ANOVAp","ANOVAq")
out2 <- out2[order(out2$ANOVAq),]
outFile <- sprintf("%s/dmr_results.txt",outDir)
write.table(out2,file=outFile,sep="\t",col=T,row=T,quote=F)

idx <- which(out2$ANOVAq < 0.2)
if (any(idx)) {
cat("-----\n")
	cat(sprintf("Loci with Q < 0.2\n"))
	print(out2[idx,,drop=F])
	cat("Locations:\n")
	sig <- rownames(out2)[idx]
	print(tgtdf[which(tgtdf$name %in% sig),])
	for (tgt in sig) {
		ph2 <- subset(ph, Var1 %in% tgt)
		p <- ggplot(ph2,aes(x=DX,y=value))
		p <- p + geom_boxplot(aes(fill=DX),lwd=1.5)
		p <- p + ggtitle(sprintf("case/control %s",tgt)) 
		p <- p + theme_bw() + theme(axis.text=element_text(size=16))
		pdf(sprintf("%s/%s.pdf",outDir,tgt))
		print(p); dev.off()
	}
cat("-----\n")
}


}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})
}

#runDMR("scz2_plt0.0001.EnhProm","glia")
#runDMR("NEURODEV","glia")
#runDMR("CaseControl","neuron")
runDMR("CaseControl","glia")
#runDMR("scz2_plt0.0001.CTCF","glia")

