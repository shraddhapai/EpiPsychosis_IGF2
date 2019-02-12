#' match SNP genotypes and those inferred from EPIC arrays

epiFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/preprocessing/rgSet_SNPcalls.txt"
epiDatFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata"
genoFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/SUS19399.hg19.sorted_-clean.EPICsnps.raw"

load(epiDatFile) #MSet.genome

geno <- read.table(genoFile,sep=" ",h=T,as.is=T)
epi <- read.delim(epiFile,sep="\t",h=T,as.is=T)

require(minfi)
require(reshape2)

# now do same for geno samples
geno <- geno[-grep("-R[23456]",geno[,2]),]
geno[,2] <- sub("-R1","",geno[,2])
geno <- geno[-which(geno[,2]=="90"),]

# now combine
epi <- cbind(epi,SNP=rownames(epi))
epim <- melt(epi); epim[,2] <- sub("X","",epim[,2])
colnames(epim)[2] <- "array_name"

geno <- geno[,c(2,7:ncol(geno))]
genom <- melt(geno) ;
colnames(genom) <- c("Sample_ID","SNP","geno")
upos <- regexpr("_",genom$SNP)
genom$SNP <- substr(genom$SNP,1,upos-1)

pd <- pData(MSet.genome)
pd <- pd[-grep("R[2345]",pd$Sample_Name),]
pd$Sample_Name <- sub("-R1","",pd$Sample_Name)
pd <- cbind(rownames(pd), pd$Sample_Name)
colnames(pd) <- c("array_name","Sample_ID")

epi2 <- merge(x=epim,y=pd,by="array_name")
both <- merge(x=genom,y=epi2,by=c("Sample_ID","SNP"))

dt <- format(Sys.Date(),"%y%m%d")
pdf(sprintf("snp2geno_%s.pdf",dt),width=6,height=8)
par(mfrow=c(5,3),mar=c(1,1,1,1))
crset <- c()
for (k in unique(both$SNP)) {
	tmp <- subset(both, SNP %in% k)
	newx <-tmp$geno+rnorm(nrow(tmp),0,0.06)
	newy <- tmp$value+rnorm(nrow(tmp),0,0.06)
	crset <- cor(newx,newy)
	plot(newx,newy,type='n',cex=0.4,
		main=sprintf("%s: %i samples",k,nrow(tmp)),
		xaxt='n',yaxt='n',cex.main=1.2)
	text(newx,newy,label=tmp$Sample_ID,cex=0.8,col=1:nrow(tmp))
}
dev.off()

