#' annotate probes by genes
#'
#' @param probe_GR (GRanges) output of getLocation(MethylSet-Genomic)
#' @param gene_GR (GRanges) gene transcript region. must have name,name2
#' 
#' @return (list of GRanges)
#' 1) promoters: probe-to-promoter matches
#' 2) genes:	probe-to-genebody matches
#' @export
annotateGene <- function(probe_GR, gene_GR, promoter_up=2000L,
	promoter_down=500L) {
	
	# promoters
	prom_GR <- promoters(gene_GR,upstream=promoter_up,
		 downstream=promoter_down)
	ol 		<- findOverlaps(probe_GR,prom_GR)
	ol_mat	<- cbind(ol@queryHits,ol@subjectHits)
	ol <- ol_mat
	anno_prom <- probe_GR[ol[,1]]
	anno_prom$promoter_name2 <- prom_GR$name2[ol[,2]]
	anno_prom$promoter_name <- prom_GR$name[ol[,2]]
	rm(ol)

	# gene bodies
	ol		<- findOverlaps(probe_GR,gene_GR)
	ol_mat	<- cbind(ol@queryHits,ol@subjectHits)
	ol 		<- ol_mat
	anno_gene	<- probe_GR[ol[,1]]
	anno_gene$gene_name2 <- gene_GR$name2[ol[,2]]
	anno_gene$gene_name <- gene_GR$name[ol[,2]]

	return(list(promoters=anno_prom,gene_body=anno_gene))
}
