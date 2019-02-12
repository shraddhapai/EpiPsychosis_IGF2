#' correlate IGF2 methylation with those of other significant loci
#' identify significant correlations
rm(list=ls())
require(minfi)

rootDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398"
mFile <- sprintf("%s/preprocessing/CaseControlEPIC_CLEAN_161221.Rdata",
				 rootDir)
dmpRes <- sprintf("%s/dmp_DX_QTL/agesexslide_PC12/dmp..dmp_DXcase.topTable.161221.FDR0.30.locations.txt",rootDir)

cat("Loading methylation data\n")
load(mFile)

cat("Loading DMP results\n")
dmp <- read.delim(dmpRes,sep="\t",h=T,as.is=T)
igf2 <- c("cg07096953","cg26401390","cg02613624")

Qthresh <- 0.3 # FDR threshold for probes to correlate

# -----------------------------------------------------------------

dmp <- subset(dmp, qval < Qthresh)

cmat <- list()
mvals <- getBeta(MSet.genome[dmp$NAME,])

for (cur in igf2) {
	cat(sprintf("%s\n",cur))
	b <- mvals[cur,]
	curcor <- sapply(1:nrow(mvals),function(x) {
				   y <- cor.test(b,mvals[x,],method="sp"); 
				   c(y$estimate,y$p.value)
	})
	rownames(curcor) <- c("cor","p.value")
	cmat[[cur]] <- t(curcor)
	rownames(cmat[[cur]]) <- rownames(mvals)
	colnames(cmat[[cur]])[1:2] <- paste(cur,colnames(cmat[[cur]])[1:2],sep=".")
}
cmat[[1]] <- cbind(cmat[[1]],qval=p.adjust(cmat[[1]][,2],method="BH"))
b <- cbind(cmat[[1]],cmat[[2]])
z <- cbind(b,cmat[[3]])
z <- cbind(z,name=rownames(z))
z <- as.data.frame(z)

locs <- getLocations(MSet.genome[rownames(z),])

# annotate with nearest gene
annoFiles 	<- read.delim("Annotations.txt",sep="\t",h=T,as.is=T)
genes <- read.delim(annoFiles[which(annoFiles[,1]=="Genes_Symbols"),3],
	sep="\t",h=F,as.is=T)
gene_GR <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),name=genes[,4],
	strand=genes[,6])

d <- distanceToNearest(locs,gene_GR)
locs$nearestGene <- gene_GR$name[subjectHits(d)]
locs$d2nearestGene <- elementMetadata(d)@listData$distance

locs_df <- as.data.frame(locs)
locs_df$name <- rownames(locs_df)
df <- merge(x=z,y=locs_df,by="name")
colnames(df)[1] <- "NAME"

dmp <- dmp[,c("NAME","pval","qval")]
colnames(dmp)[2:3] <- c("dmp.pval","dmp.qval")
df <- merge(x=df,y=dmp,by="NAME")
df <- df[order(df$cg07096953.p.value),]

options(scipen=10)
write.table(df,file="IGF2_correlations.txt",sep="\t",col=T,row=F,quote=F)

