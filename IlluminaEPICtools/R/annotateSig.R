#' annotate foreground background probes
#'
#' @details Runs variety of analyses to identify characteristics of 
#' significant probes
#' @param mset (GenomicRatioSet) 
#' @param dmpIn (DataFrame) output of minfi::dmpFinder()
#' @param annoFiles (data.frame) bed files for various annotations, 
#' with names (Type), and path to bed file (File)
#' Should have special entries for "Exons","Introns","Genes"
#' @param FisherBed (list) bed files containing genomic regions in
#' which overrepresentation of probes must be tested. If NULL, test is 
#' not performed
#' @param qCutoff (numeric 0-1) Q value cutoff for foreground
#' @param geneDomain (numeric) max distance (kb) from TSS at which a probe is considered
#' within the gene's domain. Used for GREAT analysis. Should be large enough to include
#' any probe that would be part of a gene's domain by GREAT's calculation. Should be small
#' enough that the resulting bed files can be uploaded to the GREAT server.
#' (http://bejerano.stanford.edu/great/public/html/)
#' @param outDir (char) path to output dir
#' @param fPrefix (char) a name for the analysis, to add to output files
#' e.g. fPrefix="neuron" would create output files named 
#' "dmp_neuron_annotateSig_yymmdd.pdf"
#' @param ... params for dmp_pathwayGREAT()
#' @export
annotateSig <- function(mset, dmpIn, annoFiles,FisherBed=NULL,
		qCutoff=0.05,outDir,fPrefix="",...) {
	dt <- format(Sys.Date(),"%y%m%d")
	cat(sprintf("%s: %i samples\n", fPrefix,length(mset)))
	
	fg <- rownames(dmpIn)[which(dmpIn$qval<qCutoff)]
	if (length(fg)<1) {
		cat("No significant probes.\n")
		return(TRUE)
	}
	bg <- rownames(dmpIn)
	fgpval <- rownames(dmpIn)[which(dmpIn$pval<qCutoff)]
	# extract foreground, background, get locations, make to GR.
	loc <- getLocations(mset)
	loc_fg <- loc[which(names(loc)%in% fg)]
	loc_bg <- loc[which(names(loc)%in% bg)]
	cat(sprintf("%i foreground (%i p < cutoff), %i background\n", 
				length(loc_fg),length(fgpval),length(loc_bg)))
	# force seq levels to be the same
	if (all.equal(seqlevels(loc_fg),seqlevels(loc_bg))!=TRUE) {
		cat("seqlevels not the same\n"); browser()
	# take union of both
	}

	# ---------------------------
	# Location pie charts
	# LOCATION: breakdown by chromosome
	tmp <- table(as.character(seqnames(loc_fg)))
	tmp <- tmp[paste("chr",c(1:22,c("X","Y")),sep="")]
	pdf(sprintf("%s/dmp_%s_annotateSig_LocationChrom_%s.pdf",
				outDir,fPrefix,dt),width=15,height=6)
	tryCatch({
		barplot(tmp,col=c(rep("gray70",22),"pink","cyan"),las=3,
				main=sprintf("Location of foreground probes (Q< %1.2f; N=%i)",
					qCutoff,length(loc_fg)));
	},error=function(ex){print(ex)},finally={dev.off()})

	pdf(sprintf("%s/dmp_%sannotateSig_Locations_%s.pdf",outDir,fPrefix,dt),
		width=8,height=8)
	par(cex=1.3,font.axis=1.5,las=1,lwd=3)
	tryCatch({
	# LOCATION: gene pie chart
	getLocationPie(loc_fg,annoFiles[which(annoFiles$Type %in% 
		c("Exons","Introns")),])
	title(sprintf("DMP: Q<%1.2f; N=%i; Coding/noncoding distribution",
		qCutoff,length(loc_fg)))

	# LOCATION: CGI
	cat("* Location: CGI\n")
	isl <- getIslandStatus(mset)
	isl_fg	<- isl[which(names(mset) %in% fg)]
	isl		<- table(isl_fg)
	tmp <- names(isl); isl <- as.integer(isl); names(isl) <- tmp
	pct <- round((isl/sum(isl))*100)
	names(isl) <- sprintf("%s (%s%%)", names(isl),pct)
	pal <- RColorBrewer::brewer.pal(n=4,name="YlGn")
	pie(isl,col=pal)
	title("Sig probes: Location relative to CpG density")

	# LOCATION: chromatin state
	cat("* Location: Chromatin state\n")
	cstates <- annoFiles[grep("CSTATE_",annoFiles[,1]),]
	cstates[,1] <- sub("CSTATE_","",cstates[,1])
	getLocationPie(loc_fg, cstates)
	title("Sig probes: Location in chromatin states\n(adult frontal lobe)")
	},error=function(ex) {
		print(ex)
	},finally={
		dev.off()
	})

	# annotate overlapping genes
	idx <-which(annoFiles[,1]%in%"Genes_Symbols")
	genes 	<- read.delim(annoFiles[idx,3],sep="\t",h=F,as.is=T)
	gene_GR <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),
					   name=genes[,4],strand=genes[,6])
	# extend to upstream 5K of gene
	gene_GR <- resize(gene_GR,fix="end",width=width(gene_GR)+5000)
	strand(gene_GR) <- "*"

	ol_genes<- findOverlaps(loc_fg,gene_GR)
	tmp <- data.frame(id=ol_genes@from,gene=gene_GR$name[ol_genes@to])
	tmp <- tmp[!duplicated(tmp),]
	p2g <- aggregate(tmp$gene,by=list(id=tmp$id),FUN=paste,collapse=";")
	loc_fg$gene <- rep("",length(loc_fg))
	loc_fg$gene[p2g[,1]] <- p2g[,2]

	# ---------------------------
	# FET overrepresentation
	if (!is.null(FisherBed)){
	cat("* Looking for overrepresentation\n")
	pdf(sprintf("%s/dmp_%sFET_%s.pdf",outDir,fPrefix,dt),width=15,height=5)
	tryCatch({
		tmp <- dmp_Fisher(loc_fg,loc_bg,FisherBed)
		outF <- sprintf("%s/dmp_%sFET_stats_%s.txt",outDir,fPrefix,dt)
		write.table(tmp,file=outF,sep="\t",col=TRUE,row=TRUE,quote=F)
	},error=function(ex){
		print(ex)
	},finally={
		dev.off()
	})
	}

	# ---------------------------
	# Pathway analysis input
###	cat("* Writing probes for pathway analysis\n")
###	dmp_pathwayGREAT(mset,dmpIn,geneFile=annoFiles[which(annoFiles[,1]%in% "Genes"),3],
###					 outDir=outDir,...)
###	cat("* Annotation complete\n")
	return(loc_fg)
}


