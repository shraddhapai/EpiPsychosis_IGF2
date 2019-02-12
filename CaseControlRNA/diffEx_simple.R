#' edgeR CiBERSORT diffEx Lee data

rm(list=ls())
require(edgeR)

inDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlRNAseq/SP_Lee_compare"
outDir <- sprintf("%s/STAR_SP_output",inDir)
dt <- format(Sys.Date(),"%y%m%d")

geneFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/gencode.v26.annotation.gtf.geneids_chroms.txt"
methFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlRNAseq/IGF2_meth_status_170718.txt"
idFile <- "/Users/shraddhapai/Google Drive/Projects/NARSAD_EpiPsychosis/samples/TissueSamples_updated_160922_LATEST_v2.txt"
IGF2methFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/input_files/Labrie_Brain_RNA Data-sequencing.txt"

countFile <- sprintf("%s/SCZ_RNAseq_STAR_counts.csv",inDir)
# Input was Lee's cibersort STAR data but with ensembl ids converted to
# gene symbols
# Gene signature was the Darmanis derived sig (from the autism paper)
# with the G2C synaptic genes removed.
cibOutFile <- sprintf("%s/SP_output/CIBERSORT.Output_Job9.txt",inDir)
covarFile <- sprintf("%s/covar.tsv",inDir)

library(readr)
covar <- read_delim(covarFile,"\t",escape_double=FALSE,trim_ws=TRUE)

# ---

if (!file.exists(outDir)) dir.create(outDir)

# IGF2 methylation status
igf2stat <- read.delim(IGF2methFile,sep="\t",h=T,as.is=T)
igf2stat <- igf2stat[,c("New.Code","Methylation")]
colnames(igf2stat)[1] <- "SAMPLE"
igf2stat$SAMPLE <- as.character(igf2stat$SAMPLE)
covar$SAMPLE <- sub("SC2","",covar$SAMPLE)

midx <- match(covar$SAMPLE,igf2stat$SAMPLE)
if (all.equal(igf2stat$SAMPLE[midx],covar$SAMPLE)!=TRUE) {
	cat("igf2/covar don't match"); browser()
}
covar$igf2meth <- igf2stat$Methylation[midx]

counts <- read.delim(countFile,sep=",",h=T,as.is=T)
mat_genes <- counts[,1]
counts <- counts[,-1]

#- DGEList data class
y <- DGEList(counts, group = covar$GROUP)
#- Filtering
keep <- rowSums(cpm(y)>1) >= 34
y <- y[keep, , keep.lib.sizes=FALSE]
mat_genes <- mat_genes[keep]
cat(sprintf("low count filter: %i genes survive\n", nrow(y)))
#- Normalization
y <- calcNormFactors(y,method="TMM")

cibersort_output <- read.delim(cibOutFile,sep="\t",h=T,as.is=T)
covar_ciber <- cbind(covar,cibersort_output)
# -------
pheno <- covar_ciber 
#idx <- which(pheno$neurons < 0.45)
#pheno <- pheno[-idx,]
#y <- y[,-idx]

geneFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/gencode.v26.annotation.gtf.geneids_chroms.txt"
genedef<- read.delim(geneFile,sep="\t",h=F,as.is=T)
dpos <- regexpr("\\.",genedef[,5])
genedef[,5] <- substr(genedef[,5],1,dpos-1)
midx <- match(mat_genes,genedef[,5])
genesymb <- genedef[midx,5:6]

if (all.equal(rownames(y$samples),pheno[,1])!=TRUE) {
	cat("pheno/y don't match"); browser()
}

pheno$GROUP <- factor(pheno$GROUP,levels=c("CTRL","DISEASE"))
design1 <- model.matrix(~ 1+ GROUP + SEX + AGE + PMI + neurons, data=pheno)
#- Estimate Dispersion
y1 <- estimateDisp(y,design1)
#- generalized linear model
fit1 <- glmFit(y1,design1)
#- liklihood ratio test
lrt1 <- glmLRT(fit1,coef=2)
lrt1 <- lrt1$table
lrt1$FDR <- p.adjust(lrt1$PValue,method="BH")
lrt1 <- cbind(lrt1,genesymb)
lrt1 <- lrt1[order(lrt1$FDR),]

cat(sprintf("%i genes survive FDR correction\n", sum(lrt1$FDR < 0.06)))
print(lrt1[which(lrt1$FDR < 0.06),])

write.table(lrt1,
	file=sprintf("%s/STAR_edgeR_GroupAgeSexPMINeurons_%s.txt",
	outDir,dt),sep="\t",col=T,row=F,quote=F)
require(dataExplore)
pdf(sprintf("%s/diffEx_plots_%s.pdf",outDir,dt),width=8,height=8)
plot(hist(lrt1$PValue,n=100))
drawQQplot(lrt1$PValue,main="STAR, group+sex+age+pmi+neurons")
dev.off()

# ------------------------------------------------------------
# write GSEA input
lrt1$score <- -log10(lrt1$PValue)*sign(lrt1$logFC)
tmp <- lrt1[,c("V6","score")]
tmp <- tmp[order(tmp[,2]),]
rnkFile <-sprintf("%s/STAR_edgeR_GroupAgeSexPMINeurons_%s.rnk",
	outDir,dt)
write.table(tmp,file=rnkFile,sep="\t",col=F,row=F,quote=F)

gsea <- "/Users/shraddhapai/software/gsea-3.0.jar"
gmtFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/Human_GOBP_AllPathways_no_GO_iea_October_01_2017_symbol.gmt"

gsea_cmd <- paste("java -cp ", gsea, " -Xmx1G xtools.gsea.GseaPreranked -gmx ", gmtFile, " -norm meandiv -nperm 1000 -rnk ", rnkFile, " -scoring_scheme classic -rpt_label ccRNAseq -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed 12345 -set_max 200 -set_min 10 -zip_report false -out ", outDir, " -gui false", sep="")
cat("Running GSEA\n")
print(gsea_cmd)
system(gsea_cmd) # run GSEApreranked

cat("Read GSEA results\n")
gseaDir <- sprintf("%s/%s",outDir,
	dir(path=outDir,pattern="ccRNAseq")[1])
posFile <- dir(gseaDir,pattern="gsea_report_for_na_pos")
posFile <- posFile[grep("xls",posFile)]
posFile <- sprintf("%s/%s",gseaDir,posFile)

negFile <- dir(gseaDir,pattern="gsea_report_for_na_neg")
negFile <- negFile[grep("xls",negFile)]
negFile <- sprintf("%s/%s",gseaDir,negFile)

upreg <- read.delim(posFile,sep="\t",h=T,as.is=T)
downreg <- read.delim(negFile,sep="\t",h=T,as.is=T)


#### --------------------------
#### plot IGF2 RNA by methylation
###.plotByIGF2M <- function(gn) {
###	x <- reshape2::melt(cpms[which(gene_names==gn)[1],,drop=F])
###	colnames(x)[2] <- "ID"
###	x <- x[,-1]
###	x <- merge(x=x,y=covar,by="ID")
###	x$GROUP <- factor(x$GROUP, levels=c("DISEASE","CTRL"))
###	x$DISEASE <- factor(x$DISEASE,levels=c("CTRL","SCZ","BIPOL"))
###	x$igf2meth <- factor(x$igf2meth,levels=c("l","h"))
###	p <- ggplot(x,aes(x=GROUP,y=value))
###	p <- p+geom_boxplot(aes(colour=igf2meth))+ggtitle(gn)
###	p <- p + geom_point(position=position_dodge(width=0.75),
###			aes(colour=igf2meth,group=igf2meth),cex=0.5)
###	p <- p+ylab("log-CPM")
###	print(table(x[,c("DISEASE","igf2meth")]))
###	# explain IGF2 transcription as function of IGF2 methylation, group,
###	# interaction term
###	fit <- lm(value~1+igf2meth+GROUP+igf2meth*GROUP,data=x)
###	print(summary(fit))
###	print(p)
###}
###cpms <- cpm(y,log=TRUE)
###gene_names <- lrt1$V6
###.plotByIGF2M("IGF2")

# run gsea



