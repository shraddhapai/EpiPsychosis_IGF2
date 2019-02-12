#' plot view of pvalue signal in IGF2 region
rm(list=ls())
require(GenomicRanges)
require(Gviz)
require(minfi)

# Run 
dataDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/"
dmpFile <- sprintf("%s/dmp_DX_QTL/agesexPMI_PC12_171129/dmp..dmp_DXcase.topTable.171129.txt",
	dataDir)
mData <- sprintf("%s/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata",dataDir)
annoFile <- "Annotations.txt"
trackSaveDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/geneModels"
# TF Motifs from MotifMap
motifFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/MotifMap/MotifMap_IGF2_1000_10000_161222.txt"

seqCapFile <- sprintf("/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCapEpi/methylation/IGF2_minCvg%i_tile%i.hg19coords.Rdata",minCvg,winSize)

# -----------------------------------------------------------------
# setup 

cat("* Loading methylation data\n")
load(mData)
locs <- getLocations(MSet.genome) 
cat("* Reading DMP pvalues\n")
dat <- read.delim(dmpFile,sep="\t",h=T,as.is=T)
cat("* Adding pvalues to DMP_GR")
locs <- locs[which(names(locs) %in% rownames(dat))]

midx <- match(names(locs), rownames(dat))
if (all.equal(rownames(dat)[midx],names(locs))!=TRUE) {
	cat("don't see probes\n")
}

dat <- dat[midx,]
locs$pval <- dat$P.Value
locs$qval <- dat$adj.P.Val

roi_GR <- GRanges(c("chr11","chr19"),
				 IRanges(c(2140000,15570000), c(2170000,15595000)),
				 name=c("IGF2","PGLYRP2"))
roi_GR <- roi_GR[1]
dmp_GR <- locs

#' @param dmp_GR (GRanges) probe coordinates and data to plot. First
#' metadata column will contain pvalues.
#' @param roi_GR (GRanges) windows of interest to plot. One plot will be
#' generated per such region
#' @param annoFile (char) path to Annotations.txt
	
cat("* Setting up annotations\n")
anno <- read.delim(annoFile,h=T,as.is=T,sep="\t")

genes <- read.delim(anno[grep("Genes_Symbols",anno),3],sep="\t",
					h=TRUE,as.is=TRUE)
gene_GR <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),
				   name=genes[,4],strand=genes[,6])
cat("\tGene models\n")
modelFile <- anno[grep("hg19_GeneModels",anno[,1]),3]
load(modelFile)
cat("\tDLPFC enhancers\n")
enh <- read.delim(anno[grep("CSTATE_Enh",anno[,1]),3],sep="\t",
				   h=F,as.is=T)
prom <- read.delim(anno[grep("CSTATE_Prom",anno[,1]),3],sep="\t",
				   h=F,as.is=T)
prom <- rbind(enh,prom)
prom <- subset(prom, prom[,1] %in% seqnames(roi_GR))
prom_GR <- GRanges(prom[,1],IRanges(prom[,2],prom[,3]))

cat("\tTF motifs\n")
tf <- read.delim(motifFile,sep="\t",h=T,as.is=T)
tf <- subset(tf, FDR < 0.05)
tf_GR <- GRanges(tf$chromosome, IRanges(tf$start, tf$stop),
				 strand=tf$strand,name=tf$TF.name)

cat("\tSetting up tracks\n")
itrack <- IdeogramTrack(genome="hg19")
gtrack <- GenomeAxisTrack(genome="hg19",cex=0.7)
geneTrack	<- AnnotationTrack(hg19.refseq_tx,
					genome="hg19",name="Genes",
					id=hg19.refseq_tx$geneSymbol,
					featureAnnotation="id")
displayPars(geneTrack) <- list(col="darkgreen",
				   fill="darkgreen", cex=0.7,cex.title=0.9,
				   shape="arrow",
				   arrowHeadMaxWidth=5,font=2)
promTrack <- AnnotationTrack(prom_GR,genome="hg19",name="EnhProm",
			col="orange",fill="orange",
			stacking="squish",cex.title=0.5)

cat("\tData Track\n")
dataTrack <- DataTrack(dmp_GR,data=-log10(dmp_GR$pval),type="b",
					   name="dmp p")
dmpTrack	<- AnnotationTrack(dmp_GR[which(dmp_GR$qval<0.05)],fill="red",
							   col="red",
							   name="Q<0.05",stacking="squish")
cat("\tTF track\n")
tfTrack	<- AnnotationTrack(tf_GR,fill="blue",col="red",name="TF",
							  stacking="squish",showId=TRUE,
							  id="name")

for (k in 1:length(roi_GR)) {
		chr 	<- as.character(seqnames(roi_GR)[k])
		from	<- start(roi_GR)[k]
		to		<- end(roi_GR)[k]

		tSaveFile <- sprintf("%s/%s_%s_%i_%i.Rdata",trackSaveDir,
				roi_GR$name[k],chr,from,to)
		if (!file.exists(tSaveFile)) {
				cat("Fetching from UCSC (go do some pushups)\n")
		t0 <- Sys.time()
		refGenes <- UcscTrack(genome="hg19",chromosome=chr,
							  track="refGene",
							  from=from,to=to,
							  trackType="GeneRegionTrack",
							  rstarts="exonStarts",rends="exonEnds",
							  gene="name",symbol="name2",
							  transcript="name", strand="strand",
							  fill="green",name="RefSeq Genes")
			save(refGenes,file=tSaveFile)
			print(Sys.time()-t0)
		} else {
			cat("Loading from file!\n")
			load(tSaveFile)
		
		}
		displayPars(refGenes) <- list(showId=TRUE,fill="darkgreen",
									  cex=1)

	pdfFile <- sprintf("%s/%s_detailView.pdf",dirname(dmpFile),roi_GR$name[k])
	pdf(pdfFile,width=7,height=3)
	tryCatch({
		plotTracks(list(itrack,dmpTrack,dataTrack,
						promTrack,refGenes,gtrack),
				   chromosome=chr,from=from,to=to,
				   sizes=list(0.05,0.05,0.4,0.05,0.3,0.15))
	},error=function(ex){
		print(ex)
	},finally={
		dev.off()
	})
}
###}


