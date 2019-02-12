
# gets gene expression for gene set
#' @param xpr (data.frame) master expression matrix
#' @param pd (data.frame) sample metadata
#' @param id2name (data.frame) has gene.ID->gene.name mappings
#' @param geneSet (char) names of genes to extract
#' @param aggByRegion (logical) average gene expression by all ROI
#' @return (list) 
#' 1) melted: melted data frame, grouped by age
#' 2) df: recast data.frame with rows as genes and columns as samples
getXpr <- function(xpr,pd,id2name,geneSet,aggRegions=TRUE) {
	idx <- id2name$gene.ID[which(id2name$gene.name %in% geneSet)]
	xpr <- xpr[which(rownames(xpr) %in% idx),,drop=FALSE]
	id2name <- subset(id2name, gene.ID %in% idx)

	# match row-column order - filter by samples to keep
	midx <- match(rownames(pd),colnames(xpr))
	if (all.equal(colnames(xpr)[midx],rownames(pd))!=TRUE) {
		cat("sample names don't match\n")
		browser()
	}
	xpr <- xpr[,midx,drop=FALSE] 
	xpr <- melt(xpr)
	colnames(xpr) <- c("gene.ID","geo_accession","value")

	pd <- pd[,c("geo_accession","pcw_age","regions")]
	xpr <- merge(x=xpr,y=pd,by="geo_accession") # merge with pheno
	xpr	<- merge(x=xpr,y=id2name,by="gene.ID") # merge with gene name
	xpr$gene.name <- as.character(xpr$gene.name)

	if (aggRegions) {
		agg <- aggregate(xpr$value,by=list(age=xpr$pcw_age,gene=xpr$gene.name),
			FUN=mean)
		colnames(agg) <- c("age","gene.name","value")
		df <- dcast(agg,gene.name~age,value.var="value")
	} else {
		agg <- aggregate(xpr$value,by=list(age=xpr$pcw_age,region=xpr$regions,
					gene=xpr$gene.name),FUN=mean)
		colnames(agg) <- c("age","region","gene.name","value")
		df <- dcast(agg,gene.name~age+region,value.var="value")
	}
	rownames(df) <- df[,1]; df <- df[,-1]

	return(list(melted=agg,df=df))
}

 
