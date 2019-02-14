#' mouse rnaseq analysis
rm(list=ls())

require(edgeR)
require(limma)
require(dataExplore)
require(ggplot2)
require(reshape2)

phenoFile <- "/home/shraddhapai/Epigenetics/NARSAD/input_files/MouseRNAseq2/IGF2.Mouse.RNAseq_samplekey_180625.csv"
rootDir <- "/home/shraddhapai/Epigenetics/NARSAD/input_files/MouseRNAseq2/LABV_20180620_RNA"

### run diffex based on tissue
runDiffEx <- function(gpType="all") {

ens2mgi <- read.delim("/home/shraddhapai/Epigenetics/NARSAD/anno/mm38/Ensembl_MouseGRCHm38p6_Ens2MGI_180712.txt",
	sep="\t",h=T,as.is=T)

countDir <- sprintf("%s/star",rootDir)
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/diffEx_%s/gp%s",rootDir,dt,gpType)

if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)


logFile <- sprintf("%s/diffEx.log",outDir)
sink(logFile,split=TRUE)
tryCatch({

cat('--------------------------------\n')
cat(sprintf("Group = %s\n",gpType))
cat('--------------------------------\n')

# setup
pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)

if (gpType=="FX") {
	pheno <- subset(pheno,Tissue.region=="FX")
} else if (gpType=="STR") {
	pheno <- subset(pheno,Tissue.region=="STR")
} else if (gpType=="all") {
	#pheno <- subset(pheno,Sample.ID %in% 1:20)
} else { stop("invalid gpType\n")
}

#pheno <- subset(pheno, Sample.ID %in% 1:20)

sampid <- pheno$Sample.ID

pheno <- pheno[,c("Sample.ID","Genotype","Tissue.region","Sex")]
pheno$Genotype <-factor(pheno$Genotype,levels=c("WT","KO"))
pheno$Sex <- factor(pheno$Sex,level=c("M","F"))
pheno$Tissue.region=factor(pheno$Tissue.region,level=c("FX","STR"))
cat(sprintf("Got %i samples\n",nrow(pheno)))

countFiles <- paste(sampid,"_ReadsPerGene.out.tab",sep="")
counts <- readDGE(countFiles,path=countDir,columns=c(1,4),labels=sampid,header=TRUE)
counts <- counts$counts[-1:-3,]                                                 
rownames(counts) <- gsub("\\..*","",rownames(counts))                           
write.csv(counts,file=sprintf("%s/RNAseq_STAR_counts.csv",outDir))

# make edgeR object normalize
y <- DGEList(counts, group = pheno$Genotype)
cat(sprintf("Got %i genes\n", nrow(y)))
keep <- rowSums(cpm(y)>1) >= nrow(pheno)                                                 
y <- y[keep, , keep.lib.sizes=FALSE]                                            
cat(sprintf("FILTER genes with cpm>1 in all samples = %i left\n",nrow(y)))

y <- calcNormFactors(y,method="TMM") 
cpmval <- cpm(y,normalized.lib.sizes=TRUE, log=FALSE)                              
write.csv(cpmval,file=sprintf("%s/cpm.csv",outDir))                                                   
logcpm <- cpm(y,normalized.lib.sizes=TRUE, log=TRUE)                            
write.csv(logcpm, file=sprintf("%s/logcpm.csv",outDir))

# diffex
if (gpType %in% "all") {
	cat("controlling for sex and tissue\n")
	design1 <- model.matrix(~Genotype+Sex+Tissue.region,data=pheno)
	
} else {
	cat("controlling for animal sex\n")
	design1 <- model.matrix(~Genotype+Sex,data=pheno)
}
y1 <- estimateDisp(y,design1)
fit1 <- glmFit(y1,design1)
lrt1 <- glmLRT(fit1,coef=2)
lrt1 <- lrt1$table
lrt1$FDR <- p.adjust(lrt1$PValue,method="BH")
lrt1 <- cbind(lrt1,cpmval,Gene.stable.ID=rownames(lrt1))
lrt2 <- merge(x=lrt1, y=ens2mgi,by="Gene.stable.ID")
lrt1 <- lrt2
colnames(lrt1)[which(colnames(lrt1)=="MGI.symbol")] <- "gene.name"

write.table(lrt1,file=sprintf("%s/diffEx.txt",outDir),sep="\t",col=T,row=T,quote=F)
sig <- subset(lrt1,FDR<0.05)

cat("generating diffex results\n")
png(sprintf("%s/qqplot.png",outDir))
drawQQplot(lrt1$PValue,main="Mouse RNAseq WT vs KO")
dev.off()
png(sprintf("%s/volcano1.png",outDir))
plot(lrt1$logFC,-log10(lrt1$PValue),cex=0.3,col=rgb(0,0,0,0.3),xlab="logFC",
	ylab="-log10(p)", main="Volcano, mouse WT vs KO",bty='n',las=1,cex.axis=1.3,
	cex.lab=1.2)
idx <- which(lrt1$FDR < 0.05)
points(lrt1$logFC[idx],-log10(lrt1$PValue[idx]),col='red',cex=0.3)
dev.off()
# closeup
png(sprintf("%s/volcano2.png",outDir))
plot(lrt1$logFC,-log10(lrt1$PValue),cex=0.3,col=rgb(0,0,0,0.3),xlab="logFC",
	ylab="-log10(p)", main="Volcano, mouse WT vs KO",bty='n',las=1,cex.axis=1.3,
	cex.lab=1.2,ylim=c(0,22))
idx <- which(lrt1$FDR < 0.05)
points(lrt1$logFC[idx],-log10(lrt1$PValue[idx]),col='red',cex=0.3)
dev.off()

cat(sprintf("\tDiffEx: %i with Q < 0.05\n", nrow(sig)))
sig <- sig[order(sig$PValue),]
write.table(sig,file=sprintf("%s/diffEx_Q0.05.txt",outDir),sep="\t",
	col=T,row=T,quote=F)

cat("* Running pathway analysis\n")
source("diffEx_ORA.R"); pathRes <- diffEx_ORA(lrt1)
pathRes$FET.qvalue <- p.adjust(pathRes$FET.pvalue,method="BH")
pathRes$FET.Bonf <- pathRes$FET.pvalue/nrow(pathRes)
write.table(pathRes,file=sprintf("%s/pathwayORA.txt",outDir),sep="\t",col=T,row=F,quote=F)
bonfp <- 0.05/nrow(pathRes)
sig <- subset(pathRes,FET.pvalue<bonfp)
cat(sprintf("%i pathway pass Bonferroni correction\n",nrow(sig)))
print(sig$Geneset_name)

},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})

igf2 <- "ENSMUSG00000048583"
df <- melt(logcpm) 
colnames(df)[2] <- "Sample.ID"
y <- merge(x=df,y=pheno,by="Sample.ID")

plotGene <- function(cpms, g,outFile) {
	y2 <- subset(y, Tags %in% g)
	if (nrow(y2) < 0) { cat(sprintf("%s not expressed\n", g)); return } 
	p <- ggplot(y2,aes(x=Genotype,y=value)) + geom_boxplot(aes(group=Genotype))
	p <- p + theme(text=element_text(size=30),axis.text=element_text(angle=90,size=20)) + theme_bw()
	p <- p + ggtitle(sprintf("%s: Log CPM WT vs KO",g))
	if (nrow(y2)>7) {
	z <- wilcox.test(y2$value[which(y2$Genotype %in% "WT")], y2$value[which(y2$Genotype %in% "KO")])
cat(sprintf("%s: WMW p < %1.2e\n", g,z$p.value))
}
	pdf(outFile,height=8,width=5); print(p); dev.off()
	return(p)
}

plotGene(logcpm,igf2,sprintf("%s/igf2.pdf",outDir))

if (gpType=="STR") {
	browser()
}
}

runDiffEx("FX")
runDiffEx("STR")
#runDiffEx("all")
