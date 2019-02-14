# pathway analysis mouse results
rm(list=ls())

require(org.Mm.eg.db)

gmtFile <- "/home/shraddhapai/Epigenetics/NARSAD/anno/Mouse_GOBP_AllPathways_no_GO_iea_August_01_2018_entrezgene.gmt"
gsea <- "/home/shraddhapai/software/gsea-3.0.jar"
datFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/MouseSynaptosomes/Protein_lists_for_pathway_analysis.txt"
scz_genes <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/MouseSynaptosomes/pr7b00422_si_005_valesquez et al 2017.txt"
mm2hsFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/MouseSynaptosomes/HOM_MouseHumanSequence.rpt.txt"

# prepare data
dat <- read.delim(datFile,sep="\t",h=T,as.is=T)
colnames(dat)[3] <- "protRatio"

dpos <- regexpr("\\.",dat[,1])
dat[,1] <- substr(dat[,1],1,dpos-1)

sgnChange <- rep(0, nrow(dat))
idx <- which(dat$protRatio>1); sgnChange[idx] <- 1
idx <- which(dat$protRatio<1); sgnChange[idx] <- -1

# genbank to acc numbers, which is what we have.
## Bimap interface:
x <- org.Mm.egACCNUM
# Get the entrez gene identifiers that are mapped to an ACCNUM
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
# Get the ACCNUM for the first five genes
xx[1:5]
# Get the first one
xx[[1]]
}
#For the reverse map ACCNUM2EG:
# Convert to a list
xx <- as.list(org.Mm.egACCNUM2EG)
if(length(xx) > 0){
# Gets the entrez gene identifiers for the first five Entrez Gene IDs
xx[1:5]
# Get the first one
xx[[1]]
}
acc2genbank <- cbind(names(xx),unlist(xx))
colnames(acc2genbank) <- c("Accession","genbank")

# Get MGI IDs
x <- org.Mm.egSYMBOL
# Get the entrez gene IDs that are mapped to an MGI ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
# Get the MGI gene IDs for the first five genes
xx[1:5]
# Get the first one
xx[[1]]
}
gb2symb <- cbind(names(xx), unlist(xx))
colnames(gb2symb) <- c("genbank", "symbol")

anno <- merge(x=gb2symb, y=acc2genbank,by="genbank")


#orthologues from jax
mm2hs <- read.delim(mm2hsFile,sep="\t",h=T,as.is=T)
mouse  <- mm2hs[seq(1,nrow(mm2hs),2),"Symbol"]
human  <- mm2hs[seq(2,nrow(mm2hs),2),"Symbol"]
mm2hs <- cbind(mouse=mouse,human=human)

# scz genes from velasquez
scz <- read.delim(scz_genes,sep="\t",h=T,as.is=T)
colnames(scz) <- scz[1,]; scz <- scz[-1,]
scz <- scz[,2]

out <- c()
for (k in scz) {
	x <- trimws(unlist(strsplit(k,"\\\\")))
	idx <- grep("Gname=",x)
	if (any(idx)) {
		out <- c(out, sub("Gname=","",x[idx]))
	}
}

y <- merge(x=dat,y=anno,by="Accession")
y <- y[,c("symbol","Significance")]
y <- y[!duplicated(y),]

cat(sprintf("Of %i, entries, %i survive\n",nrow(dat),nrow(y))) 

dat <- y
dat$Pval <- 10^(-dat$Significance)

