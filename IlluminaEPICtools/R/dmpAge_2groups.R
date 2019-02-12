#' compare <15 yo to >15yo. probe-wise analysis
#'
#' @param pd (data.frame) phenotype matrix
#' @param mset (MethylSet-Genomic) 
#' @param outDir (char) path to dir where results should be written
#' @param runORA (logical) use gprofiler to run functional ORA
#' @param Qcutoff (numeric between 0 and 1) Q value for sig probes
#' @param gene_GR (GRanges) gene definitions, used only if runORA=TRUE
#' @param runGroups (char) which groups to run for? [neuron|glia]
#' @import gProfileR
#' @import gplots
#' @return Nothing. Side effect of writing following tables/plots:
#' 1) dmp_age_2groups_<tis>_<dt>: stats for significant probes. From dmpRun()
#' 2) dmp_age_2groups_<tis>_<dt>.pdf: nominal p hist, qqplot, intercept hist.
#' from dmpRun()
#' 3) dmp_age_<tis>_heatmap_<dt>: heatmap of sig probes.
#' 4) age_2groups_<tis>_[ALL|UP|DOWN][prom|gene]_gProfileR_<dt>: gprofiler
#' output if runORA=TRUE
#' @export
dmpAge_2groups <- function(pd, mset,outDir,runORA=FALSE,
		Qcutoff=0.05,gene_GR,runGroups="neuron") {
dt <- format(Sys.Date(),"%y%m%d")

logFile <- sprintf("%s/dmpAge_Lenient_%s.log",outDir,dt)
sink(logFile,split=TRUE) 

tryCatch({

idx <- which(pd$age_group %in% 5)
pd <- pd[-idx,]
mset <- mset[,-idx]
M <- getBeta(mset)
cat(sprintf("Excluding adult sample: %i samples left\n",nrow(pd)))

pd$age_gp2		<- rep("older",nrow(pd))
pd$age_gp2[which(pd$age_group %in% c(1,2))] <- "younger"
gp <- factor(pd$age_gp2,levels=c("younger","older"))


.runGP <- function(fg,oF) {
	cat("\t\tRunning gprofiler\n")
	gpres <- gProfileR::gprofiler(fg, organism="hsapiens",
		   max_p_value=0.05,max_set_size=500,correction_method="fdr",
		   custom_bg=bg_prom)
	gpres <- gpres[order(gpres$p.value,decreasing=FALSE),]
	cat("\t\tWriting results to file\n")
	write.table(gpres,file=oF,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
	oF <- sub(".txt","_short.txt",oF)
	gpres <- gpres[,c("p.value","term.id","domain",
					  "term.name","intersection")]
	write.table(gpres,file=oF,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}

# run dmp on tissue-wise basis
gpTypes <- list(
	neuron=which(pd$tissue %in% "NeuN+"),
	glia=which(pd$tissue %in% "NeuN-")
)

for (runType in runGroups) {
	# run dmp
	outFile <- sprintf("%s/dmp_age_2groups_%s_%s",outDir,runType,dt)
	
	curr_idx 	<- gpTypes[[runType]]
	curr_pd		<- pd[curr_idx,]
	out <- IlluminaEPICtools::dmpRun(mset[,curr_idx],gp[curr_idx],
				  outFile,sprintf("%s: Young vs Old",runType),
				  Qcutoff=Qcutoff)
	#plotCpg(mset, cpg=as.character(down), pheno=gp)
	
	probe_GR <- GRanges(out$CHROM,IRanges(out$START_POS,out$END_POS),
					name=as.character(out$PROBE))

	qval <- out$qval
	if (any(qval < Qcutoff)) {
	
	sig		<- subset(out, out$qval < Qcutoff)
	curr_M	<- M[which(rownames(M) %in% rownames(sig)),curr_idx]
	if (sum(out$qval < Qcutoff) > 5) {
	pdf(sprintf("%s/dmp_age_%s_heatmap_%s.pdf",outDir,runType,dt))
	heatmap.2(curr_M,trace='none',
			  labCol=curr_pd[,"age"],scale='row')
	dev.off()
	} else {
		cat("\ttoo few for heatmap\n")
	}

	if (runORA){

	annoGene <- annotateGene(probe_GR,gene_GR)

	prom		<- annoGene[["promoters"]]
	bg_prom		<- unique(prom$promoter_name2) #background

	gene		<- annoGene[["gene_body"]]
	bg_gene		<- unique(gene$gene_name2)

	cat("* Running gProfileR for sig results...\n")
	cat("**** Promoters ****\n")
	# all sig
	sig_idx <- which(as.character(prom$name) %in% sig$PROBE)
	sig_gene<- unique(prom$promoter_name2[sig_idx])
	cat(sprintf("\tAll sig: N=%i probes => %i promoters\n",
				length(sig_idx),length(sig_gene)))
	.runGP(sig_gene, sprintf("%s/age_2groups_%s_ALLSIGprom_gProfileR_%s.txt",
							 outDir,runType,dt))
	
	# goes down
	down<-sig$PROBE[which(sig$intercept < 0)]
	down_idx	<- which(as.character(prom$name)%in%down)
	down_gene	<- unique(prom$promoter_name2[down_idx])
	cat(sprintf("\tDown sig: N=%i probes => %i promoters\n",
				length(down_idx),length(down_gene)))
	.runGP(down_gene, sprintf("%s/age_2groups_%s_DOWNprom_gProfileR_%s.txt",
							  outDir,runType,dt))
	
	# goes up
	up<-sig$PROBE[which(sig$intercept > 0)]
	up_idx	<- which(as.character(prom$name)%in%up)
	up_gene	<- unique(prom$promoter_name2[up_idx])
	cat(sprintf("\tUp sig: N=%i probes => %i promoters\n",
				length(up_idx),length(up_gene)))
	.runGP(up_gene, sprintf("%s/age_2groups_%s_UPprom_gProfileR_%s.txt",
							outDir,runType,dt))

	cat("**** Gene body ****\n")
	# all sig
	sig_idx <- which(as.character(gene$name) %in% sig$PROBE)
	sig_gene<- unique(gene$gene_name2[sig_idx])
	cat(sprintf("\tAll sig: N=%i probes => %i genes\n",
				length(sig_idx),length(sig_gene)))
	.runGP(sig_gene, sprintf("%s/age_2groups_%s_ALLSIGgenes_gProfileR_%s.txt",
							 outDir,runType,dt))
	
	# goes down
	down<-sig$PROBE[which(sig$intercept < 0)]
	down_idx	<- which(as.character(gene$name)%in%down)
	down_gene	<- unique(gene$gene_name2[down_idx])
	cat(sprintf("\tDown sig: N=%i probes => %i genes\n",
				length(down_idx),length(down_gene)))
	.runGP(down_gene, sprintf("%s/age_2groups_%s_DOWNgenes_gProfileR_%s.txt",
							  outDir,runType,dt))
	
	# goes up
	up<-sig$PROBE[which(sig$intercept > 0)]
	up_idx	<- which(as.character(gene$name)%in%up)
	up_gene	<- unique(gene$gene_name2[up_idx])
	cat(sprintf("\tUp sig: N=%i probes => %i genes\n",
				length(up_idx),length(up_gene)))
	.runGP(up_gene, sprintf("%s/age_2groups_%s_UPgenes_gProfileR_%s.txt",
							outDir,runType,dt))
	}
	} else {
		cat(sprintf("no probes with Q < %1.2f\n",Qcutoff))
	}
	}
	
},error=function(ex){
	print(ex)
},finally={
		cat("closing log\n")
	sink(NULL)
})
}
