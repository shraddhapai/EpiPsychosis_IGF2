#' hypergeometric
rm(list=ls())
require(org.Mm.eg.db)

#' Pathway ORA for synaptosomal proteins relative to all others in proteome
#' Test to ensure that synaptosomal enrichment worked.

datFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/MouseSynaptosomes/Protein_lists_for_pathway_analysis.txt"

# prepare data
dat <- read.delim(datFile,sep="\t",h=T,as.is=T)
require(dataExplore)
pdf("peptide_qqplot.pdf"); drawQQplot(dat[,3],main="Synaptosomal proteins",
	cex.axis=1.5,bty='n')
dev.off()
dpos <- regexpr("\\.",dat[,1])
dat[,1] <- substr(dat[,1],1,dpos-1)

cat("getting entrez\n")
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
tmp <- sapply(names(xx), function(nm) {
	ln <- length(xx[[nm]])
	cbind(rep(nm,ln),xx[[nm]])
})

acc2genbank <- do.call("rbind",tmp)
colnames(acc2genbank) <- c("genbank","Accession")

cat("getting mgi\n")

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
# now map accession to symbols
y <- merge(x=dat,y=anno,by="Accession")
y <- y[,c("symbol","Significance")]
colnames(y)[1] <- "Gene.Symbol"
y <- y[!duplicated(y),]

# run fet
dat <- y
dat$Pval <- 10^(-dat$Significance)

diffEx <- dat; rm(dat)

#----------------------------------------
# pathway analysis

pathGMT <- "/home/shraddhapai/Epigenetics/NARSAD/anno/Mouse_AllPathways_September_01_2018_symbol.gmt"
allpath <- netDx::readPathways(pathGMT,MIN=15,MAX=200)
diseaseGMT <- "~/Epigenetics/NARSAD/anno/mouse/MGI_DO.gmt"
disease <- netDx::readPathways(diseaseGMT,MIN=10,MAX=200)

allpath <- c(allpath,disease)

res <- matrix(NA,nrow=length(allpath), ncol=6)
ctr <- 1
fg_genes <- c()
is_sig <- unique(diffEx$Gene.Symbol)
all_prot <- unique(unlist(allpath))
not_sig <- setdiff(all_prot,is_sig)
for (nm in names(allpath)) {
	curgroup <- allpath[[nm]]
	fg_in_group <- intersect(is_sig, curgroup)
	bg_in_group <- intersect(not_sig,curgroup)
	mat <- matrix(NA,nrow=2,ncol=2)
	mat[1,1] <- length(fg_in_group)
	mat[1,2] <- length(is_sig)-length(fg_in_group)
	mat[2,1] <- length(bg_in_group)
	mat[2,2] <- length(not_sig)- length(bg_in_group)
	x <- fisher.test(mat)
	cat(sprintf("%s: %i total: %i in sig (p< %1.2e)\n",
			nm,length(curgroup),length(fg_in_group),x$p.value))
	
	res[ctr,] <- c(length(curgroup),length(fg_in_group),length(bg_in_group),
		length(fg_in_group)/length(is_sig), length(bg_in_group)/length(not_sig),
		x$p.value)
	fg_genes <- c(fg_genes,paste(fg_in_group,collapse=","))
	
	ctr <- ctr+1
}


colnames(res) <- c("set_size","FG.ol","BG.ol",
		"%FG","%BG","FET.pvalue")
res[,4:5] <- round(res[,4:5]*100,digits=2)
res <- as.data.frame(res)
res$FG_genes <- fg_genes
res$Geneset_name <- names(allpath)
res <- subset(res, res$FG.ol>=1 & res$BG.ol>=1)
cat(sprintf("Filter path with 1+ FG and 1+ BG => %i pathways\n",nrow(res)))
res <- res[,c(ncol(res),1:(ncol(res)-1))]
cat(sprintf("Bonferroni p is: %1.2e\n", 0.05/nrow(res)))

fdr <- p.adjust(res$FET.pvalue, method="BH")
res$qval <- fdr
res <- res[order(res$qval),]
print(head(res))
write.table(res,file="synaptosome_enrichment.txt",sep="\t",col=T,row=F,quote=F)

# write gmt file
outFile <- "synaptosome_pathenrichment.gmt"
system(sprintf("cat /dev/null > %s",outFile))

res <- subset(res,qval < 0.05)
cat(sprintf("%i pathways with Q < 0.05\n",nrow(res)))

for (k in res$Geneset_name) {
	gn <- allpath[[k]]
	cat(sprintf("%s\t%s\t%s\n",k,k,paste(gn,collapse="\t")),file=outFile,
		append=TRUE)
}


