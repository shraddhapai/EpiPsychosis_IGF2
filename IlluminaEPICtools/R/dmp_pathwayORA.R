#' pathway overrepresentation by FET/binomial
#'
#' @param mset (GenomicRatioSet)
#' @param dmpIn (DataFrame) output of dmp-wise analysis.
#' @param bg_GR (GRanges) background
#' @param gene_GR (GRanges) gene bodies. Gene symbol should be in "name"
#' slot
#' @param geneDomain (integer) distance (kb) upstream from TSS
#' TSS that should be included in gene's domain. e.g.  Setting 
#' geneDomain=10 would result in [TSS-10kb,TES] being counted as a gene's
#' domain. 
#' @param pathwayFile (char) path to gmt file with gene-sets
#' @param outDir (char) path to dir to write output files
#' @param selMode (char) pvalue|qvalue; if pvalue, foreground is probes
#' with p < 0.05; if qvalue, foreground is probes with Q < 0.05
#' @param gmtFDRcutoff (numeric 0-1) Q-value cutoff to use for writing gmt
#' @param verbose (logical) print messages
#' @import GenomicRanges
#' @return (data.frame) pathway statistics including Name, num probes 
#' overlapping in foreground (n_fg), in background (n_bg), pvalues for 
#' Hypergeometric (Hypergeom_p) and Binomial two-sided test (Binom_p).
#' and Benjamini-Hochberg adjusted Q values too (Hypergeom_Q, Binom_Q)
#' @export
dmp_pathwayORA <- function(mset,thresh=0.05,dmpIn,gene_GR,geneDomain=5,
	pathwayFile,outDir,selMode="pvalue",gmtFDRcutoff=0.05,verbose=TRUE) {

	cat("Pathway analysis\n")
	cat(sprintf("\t %s < %1.4f\n", selMode, thresh))
	dt <- format(Sys.Date(),"%y%m%d")
	if (selMode %in% "pvalue")
		fg <- rownames(dmpIn)[which(dmpIn$pval<thresh)]
	else
		fg <- rownames(dmpIn)[which(dmpIn$qval<thresh)]
			
	bg <- rownames(dmpIn)
	loc <- getLocations(mset)
	fg_GR <- loc[which(names(loc)%in% fg)]
	bg_GR <- loc[which(names(loc)%in% bg)]
	cat(sprintf("%i foreground, %i background\n", 
				length(fg_GR),length(bg_GR)))

	pathwayList <- readPathways(pathwayFile)
##	warning("**** cutting pathway list to 500") # for debugging
##	pathwayList <- pathwayList[1:500]

	# define gene domain
	cat(sprintf("* Gene domain is [TSS-%ikb,TES]\n", 
				geneDomain))
	gene_GR <- resize(gene_GR,width=width(gene_GR)+(geneDomain*1000),
					  fix="end")

	lenf <- length(fg_GR)
	lenb <- length(bg_GR)
	##out <- matrix(nrow=length(pathwayList),ncol=5)
	t0 <- Sys.time()
	out <- mclapply(1:length(pathwayList), function(k) { 
		cur <- gene_GR[which(gene_GR$name %in% pathwayList[[k]])]
		if (verbose) cat(sprintf("%s\n\t%i genes: ",names(pathwayList)[k],
				length(cur)))
		fg <- subsetByOverlaps(fg_GR,cur); nf <- length(fg)
		bg <- subsetByOverlaps(bg_GR,cur); nb <- length(bg)
		if (verbose) cat(sprintf("%i fg, %i bg\n", nf,nb))

		mat <- matrix(nrow=2,ncol=2)
		mat[1,] <- c(nf, lenf-nf) # fg - yes, no
		mat[2,] <- c(nb, lenb-nb)

		ol <- findOverlaps(fg_GR,cur)
		tmp <- c(c(nf,nb),
			fisher.test(mat)$p.value, # hypergeometric
			binom.test(nf,lenf,p=nb/lenb)$p.value,
			paste(sort(unique(cur$name[ol@to])),collapse=";"))
		tmp
	},mc.cores=4)
	t1 <- Sys.time()
	print(t1-t0)

	# write result
	out <- matrix(unlist(out),byrow=TRUE,ncol=5)
	colnames(out) <- c("n_fg","n_bg","Hypergeom_p","Binom_p","genes_fg")
	out <- as.data.frame(out)
	out$Hypergeom_p <- as.numeric(as.character(out$Hypergeom_p))
	out$Binom_p <- as.numeric(as.character(out$Binom_p))
	out$Hypergeom_Q <- p.adjust(out$Hypergeom_p,method="BH")
	out$Binom_Q <- p.adjust(out$Binom_p,method="BH")
	out$Pathway <- names(pathwayList)

	outF <- sprintf("%s/dmp_pathways_stats_%s.%1.2f.%s.txt",
					outDir,selMode,thresh,dt)
	out <- out[,c("Pathway","n_fg","n_bg","Hypergeom_p","Binom_p",
				  "Hypergeom_Q","Binom_Q","genes_fg")]
	out <- out[order(out$Hypergeom_Q),]
	write.table(out,file=outF,sep="\t",col=T,row=F,quote=F)

	# write gmt of significant pathways
	idx <- which(out$Hypergeom_Q<gmtFDRcutoff)
	if (any(idx)) {
		cat(sprintf("HyperG: # pathways Q < %1.2f = %i\n",gmtFDRcutoff,
					length(idx),gmtFDRcutoff))
		outF <- sprintf("%s/dmp_pathways_hyperG_Q%1.2f_%s.gmt",outDir,
					gmtFDRcutoff,dt)
		if (file.exists(outF)) unlink(outF)
		system(sprintf("touch %s",outF))
		curp <- which(names(pathwayList)%in% out$Pathway[idx])
		ol <- findOverlaps(fg_GR,gene_GR)
		fg_genes <- unique(gene_GR$name[ol@to])
		for (nm in names(pathwayList)[curp]) {
			# keep only genes contributing to pathway
			g <- intersect(pathwayList[[nm]], fg_genes)
			cat(sprintf("%s\t%s\t%s\n",nm,nm,paste(g, collapse="\t")),
				file=outF,append=TRUE)
		}
	}

	# repeat for binom
	idx <- which(out$Binom_Q<gmtFDRcutoff)
	if (any(idx)) {
		cat(sprintf("Binomial: # pathways Q < %1.2f = %i\n",gmtFDRcutoff,
					length(idx)))
		outF <- sprintf("%s/dmp_pathways_binom_Q%1.2f_%s.gmt",outDir,
			gmtFDRcutoff,dt)
		if (file.exists(outF)) unlink(outF)
		system(sprintf("touch %s",outF))
		curp <- which(names(pathwayList)%in% out$Pathway[idx])
		ol <- findOverlaps(fg_GR,gene_GR)
		fg_genes <- unique(gene_GR$name[ol@to])
		for (nm in names(pathwayList)[curp]) {
			# keep only genes contributing to pathway
			g <- intersect(pathwayList[[nm]], fg_genes)
			cat(sprintf("%s\t%s\t%s\n",nm,nm,
				paste(g,collapse="\t")),file=outF,append=TRUE)
		}
	}

	out
}
