#' annotate diffEx genes from CaseControl RNAseq
rm(list=ls())
source("../../util/readPathwayFile.R")
source("../../util/cleanPathwayName.R")

# ---------------------------------------------------------------
outDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlRNAseq/SP_Lee_compare/STAR_SP_output"#/diffEx_CIBERSORT_incXM_171024"

diffEx <- sprintf("%s/STAR_edgeR_GroupAgeSexPMINeurons_180111.txt",outDir)

gseaResNeg <- sprintf("%s/ccRNAseq.GseaPreranked.1508871379422/gsea_report_for_na_neg_1508871379422.xls",outDir)

#bspan <- "/Users/shraddhapai/Google Drive/genome_annotation/BrainSpan/BrainSpan_TableS13_stringent.gmt"
bspan <- "/Users/shraddhapai/Google Drive/genome_annotation/BrainSpan/BrainSpan_TableS13.gmt"

gmtFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/Human_GOBP_AllPathways_no_GO_iea_October_01_2017_symbol.gmt"

# ---------------------------------------------------------------

diffEx_full <- read.delim(diffEx,sep="\t",h=T,as.is=T)

# Test ORA in a variety of genesets
source("diffEx_ORA.R")
ora_res <- diffEx_ORA(diffEx_full)
write.table(ora_res,file=sprintf("%s/diffEx_ORA.txt",outDir),sep="\t",
	col=T,row=F,quote=F)

# extract synaptic genes in +ve list
###gseaRes <- read.delim(gseaResNeg,sep="\t",h=T,as.is=T)
###gseaRes <- subset(gseaRes,FDR.q.val < 0.05)
# get diffEx genes that contribute to brain-related pathway enrichment
###synPath <- gseaRes$NAME[grep("GABA|AMPA|CHANNEL|REELIN|MUSCARINIC|KAINATE|GLUTAMATE|AXON|DENDRIT|SYNAP|NEURO|LEARN|CAMKII|ACTION|BEHAVIO",gseaRes$NAME)]
###cat(sprintf("%i of %i pathways are synapse-related\n",
###	length(synPath),nrow(gseaRes)))
###path <- readPathways(gmtFile,getOrigNames=TRUE)
###path <- path[which(names(path) %in% synPath)]
#### find the genes that contribute to this signal
###sig_genes <- diffEx_full$gene.name[which(diffEx_full$PValue<0.05)]
###gn <- unlist(path); names(gn)<- NULL
###syngenes <- unique(intersect(gn,sig_genes))
###cat(sprintf("\t%i of %i sig genes are in synaptic pathways\n", 
###	length(syngenes),length(sig_genes)))
###write.table(syngenes,file="diffExp0.05_GenesInNeuroPathways.txt",
###	sep="\t",col=F,row=F,quote=F)


# plot diffEx genes in brain development (BrainSpan)
source("../../ExtSources/BrainSpan_plotXpr.R")
source("../../util/summarySE.R")
source("../../util/multiplot.R")
bspan <- readPathways(bspan,MIN=0,MAX=10000)

bspan <- bspan[grep("Synap|Potas|Dendri|Myeli",names(bspan))]

#### Plot PC1 of geneset with BrainSpan marker sets
#### ----------------------------------------------
outBSpan <- list()
for (nm in names(bspan)[grep("Synapse_development",names(bspan))]){
	print(nm)
	sig_genes <- diffEx_full$gene.name[which(diffEx_full$PValue<0.05)]
					
	outBSpan[[nm]] <- BrainSpan_plotXpr(sig_genes,
		ttl="diffEx genes in psychosis",
		markerGenes=bspan[[nm]],
		markerLabel=nm,outDir=sprintf("%s/BrainSpan",outDir),regions2show="prefrontal")
}

compPC <- list() # compile in one list and then melt
tmp <- outBSpan[[1]]$data_PC1[,c("age","PC1")]; 
colnames(tmp)[2] <- "PC1"
tmp$set <- "data";
compPC[["data"]] <-tmp
for (k in 1:length(outBSpan)) {
	tmp <- outBSpan[[k]]$markerPC; 
	setName <- sub("Neuron_development:","",outBSpan[[k]]$markerLabel)
	tmp$set <- substr(setName,1,14)	
	rownames(tmp) <- NULL
	compPC[[outBSpan[[k]]$markerLabel]] <- tmp
}
compPC <- do.call("rbind",compPC)

p <- ggplot(compPC,aes(x=age,y=PC1,colour=set))+geom_line(lwd=2)
comp2 <- subset(compPC,set=="data")
p <- p + geom_line(data=comp2,aes(x=age,y=PC1),colour="black",lwd=2)
# --- stuff that should probably be somewhere else
p <- p+ylim(c(0,120))+coord_trans(x="log10")
p <- p+xlab("Post-conception age (days)")+ylab("PC1 (% max)")
p <- BS_addRefs(p)
p <- p + theme(axis.text=element_text(size=16))
pdf(sprintf("BrainSpan_plotPC.pdf"),width=8,height=4)
print(p)
dev.off()

# compile PC1 correlations
corBSpan <- lapply(outBSpan,function(x) x$cr)
corBSpan <- lapply(corBSpan,function(x) { cbind(x$estimate,x$p.value)})
corBSpan <- do.call("rbind",corBSpan); rownames(corBSpan) <- NULL
corBSpan <- data.frame(corBSpan)
colnames(corBSpan) <- c("cor","cor.p")
corBSpan$BSpan_set <- names(bspan)
write.table(corBSpan[,c(3,1,2)],file=sprintf("%s/BrainSpan/correlations.txt",
	outDir),sep="\t",col=T,row=F,quote=F)

plotList <- lapply(outBSpan,function(x) x$plot)
pdf(sprintf("%s/BrainSpan/markerSet_correlation.pdf",outDir),
		width=11,height=11)
tryCatch({
	nr=3;nc=2;
  for (k in seq(1,length(plotList),nr*nc)) {
    print(multiplot(plotlist=plotList[k:(k+((nr*nc)-1))],
        layout=matrix(1:(nr*nc),nrow=nr,ncol=nc)))
  }
},error=function(ex){print(ex)},finally={dev.off()})

#### BrainSpan - plot genes of interest
BrainSpan_plotXpr(c("IGF2"), outDir=outDir) #,"PGLYRP2","SYN1","F13A1","COL6A3","TH"),
#outDir=outDir)
								



