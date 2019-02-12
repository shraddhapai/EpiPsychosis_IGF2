#' given a set of genes, plot their expression through
#' human brain development using the BrainSpan dataset.
#' @param geneNames (char) vector of gene names to plot
#' @param ttl (char) title for plot
#' @param regions2show (char) one of "all_frontal","prefrontal"; controls
#' which brain regions are included in extracted expression
#' @param markerGenes (char) vector of marker genes to add. In the PC
#' plot, PC will separately be perofr
#' @param outDir (char) dir to save plots
#' @param numPerm (integer) adds non-parametric pvalue if not NA
#' @return plots
BrainSpan_plotXpr <- function(geneNames,ttl="gene set",
	regions2show="prefrontal",markerGenes=NULL,markerLabel,outDir=".",
	numPerm=NA){

regionSets <- list(
	all_frontal=c("NCX","OFC","DFC","VFC","MFC","M1C","MGE","LGE","CGE"),
	prefrontal=c("DFC"),# "CGE","LGE","MGE"),
	mix=c("NCX","HIP","CBC")
)

require(reshape2)
require(Biobase)
require(ggplot2)

dt <- format(Sys.Date(),"%y%m%d")

### Prepare expression data - run once.
###bsFile <- "/Users/shraddhapai/Google Drive/genome_annotation/BrainSpan/GSE25219-GPL5175-ExprSet.Rdata"
###cat("* Loading BrainSpan\n")
###load(bsFile)
###exprSet <- GSE25219_GPL5175_ExprSet; rm(GSE25219_GPL5175_ExprSet)
###
###cat("* BS: Preparing sample data\n")
###pd <- pData(exprSet)[,1:14]
###xpr <- exprs(exprSet)
###fd <- featureData(exprSet)@data
###
###cat("* BS: Extracting gene names (takes 10 sec)\n")
###gn <- strsplit(as.character(fd$gene_assignment),split="///")
###gn2 <- list()
###for (k in 1:length(gn)) {
###	blah <- strsplit(gn[[k]],"//")
###	gn2[[k]] <- unique(unlist(lapply(blah,function(x) x[2])))
###}
###gn2 <- lapply(gn2,trimws)
###names(gn2) <- as.character(fd$ID)
###gn2 <- melt(gn2)
###colnames(gn2) <- c("gene.name","ID")
###gn2$ID <- as.integer(gn2$ID)
###
###geneID2name <- gn2
bsFile <- "~/Epigenetics/NARSAD/anno/BrainSpan_prepared.Rdata"
###save(geneID2name, pd, xpr, fd, file=bsFile)

load(bsFile)
gn2 <- geneID2name
gn2$gene.name <- as.character(gn2$gene.name)
gn2$gene.ID <- as.integer(gn2$gene.ID)
geneID2name <- gn2

# Recode age into PCW-days 
pd$pcw_age <- BS_convertAge(pd$characteristics_ch1.4)
# filter by regions to keep
pd$regions <- as.character(pd$characteristics_ch1.1)
pd$regions <- sub("region: ", "",pd$regions)
tokeep <- which(pd$regions %in% regionSets[[regions2show]])
pd <- pd[tokeep,]

# -------------------------------------------------------
# Plot heatmap & PC1 of geneset
# -------------------------------------------------------
# heat map of regions
if (length(geneNames)>=6) {
	# get xpr of genes of interest
	geneSet <- getXpr(xpr,pd,geneID2name,geneNames,aggRegions=TRUE)
cat("y")
	agg2 <- geneSet[[2]]
	pr <- prcomp(agg2,scale=TRUE)
	pc1 <- data.frame(age=as.integer(colnames(agg2)),PC1=pr$rotation[,1],
			PC2=pr$rotation[,2])

	# do same with marker genes
cat("going for markers\n")
	markers <- getXpr(xpr,pd,geneID2name,markerGenes,aggRegions=TRUE)	
	aggmark2 <- markers[[2]]
	if (nrow(aggmark2)>=2) {
	prmark <- prcomp(aggmark2,scale=TRUE)
	pcmark <-data.frame(age=as.integer(colnames(aggmark2)),
		PC1=prmark$rotation[,1]) 
	cr <- cor.test(pc1$PC1,pcmark$PC1)
	} else {
	pcmark <-data.frame(age=as.integer(colnames(aggmark2)),
		PC1=as.numeric(aggmark2[1,]))
	cr <- cor.test(pc1$PC1,as.numeric(aggmark2[1,]))
	}

	if (!is.na(numPerm)) {
		cat("computing correlation for non-parametric tests\n")
		set.seed(123); # make reproducible
		require(doParallel)
		cl <- makeCluster(12)
		registerDoParallel(cl)
		t0 <- Sys.time()
		perm_cor <- foreach (perm_ctr=1:numPerm,.packages="reshape2") %dopar% {
				samp <- sample(geneID2name[,1],length(geneNames),TRUE)
				source("../../ExtSources/BrainSpan_getXpr.R")
				
				tmp <- getXpr(xpr,pd,geneID2name,samp,aggRegions=TRUE)
				prtmp <- prcomp(tmp[[2]],scale=TRUE)
				# get correlation of randomly-sampled genes with marker genes
			  return(cor.test(prtmp$rotation[,1],pcmark$PC1)$estimate)
		}
		t1 <- Sys.time()
		cat(sprintf("Took %1.2f seconds\n", t1-t0))
		perm_cor <- unlist(perm_cor)
		cat("*********\n")
		permp <- sum(perm_cor> cr$estimate)/length(perm_cor)
		if (permp < 1/length(perm_cor)) permp <- 1/length(perm_cor)
		cat(sprintf("Correlation p-value from random-sample test:\n"))
		cat(sprintf("Real cor: %1.2f ; median null cor: %1.2f; (%i of %i), p < %1.2f\n",
				cr$estimate,median(perm_cor), 
				sum(perm_cor > cr$estimate), length(perm_cor),permp))
		cat("distr. of null cor\n")
		print(summary(perm_cor))
		cat("*********\n")
			stopCluster(cl)
	}

	

	.rescale <- function(x) { 
		x <- (x-min(x,na.rm=T))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
		x*100
	}

	pc1$PC1 <- .rescale(pc1$PC1)
	pcmark$PC1 <- .rescale(pcmark$PC1)

	p <- ggplot(pc1,aes(x=age,y=PC1)) + geom_point() 
	p <- p + geom_line(colour="black",lwd=2)
	#	geom_smooth(method="loess",span=0.1,se=FALSE,colour="black",lwd=2)
	p <- p+xlab("Post-conception age (days)")+ylab("PC1 of gene-set")
	p <- p+ggtitle(sprintf("%s:PC1:N=%i genes (cor=%1.2f; p < %1.2e) {%s}",
		ttl, length(unique(geneSet[[1]]$gene.name)),cr$estimate,cr$p.value,
		paste(regionSets[[regions2show]],collapse=",")))


	# add marker line
	p <- p + geom_line(data=pcmark,aes(x=age,y=PC1),colour="red")
	#p <- p + #geom_smooth(data=pcmark,aes(x=age,y=PC1),method="loess",
			#colour="red",span=0.1,se=FALSE)
	p <- p + annotate("text",label=markerLabel,
			x=70*365,y=115,colour="red",hjust=1,fontface=3)
	# add dev milestone references
	p <- BS_addRefs(p)
	p <- p+ylim(c(0,120))+coord_trans(x="log10")
	
	pdf(sprintf("%s/BrainSpan_PCplot_%s_%s.pdf",outDir,markerLabel,dt),
		width=8,height=2.5)
	tryCatch({
		p <- p + theme(axis.text=element_text(size=16))
		print(p)
	},error=function(ex){print(ex)},finally={dev.off()})

	if (!is.null(markerGenes)) {
				return(list(cr=cr,plot=p,
							data_PC1=pc1,markerPC=pcmark,
							markerLabel=markerLabel))
   }
}

# -------------------------------------------------------
# plot individual gene trends
# -------------------------------------------------------
if (length(geneNames)<=10) {
### Next step : Plot mean+sem per gene through dev
### Then do a heatmap of the genes, and colour code by dev age 
	geneSet <- getXpr(xpr,pd,geneID2name,geneNames,aggRegions=FALSE)
	sel_xpr <- geneSet[[1]]
plotList <- list()
for (curgene in unique(sel_xpr$gene.name)){
	print(curgene)
	cur <- subset(sel_xpr, gene.name %in% curgene)
	x <- summarySE(cur,measurevar="value",groupvars=c("age","region"))
	p <- ggplot(x, aes(x=age, y=value,colour=region))
	p <- p +geom_errorbar(aes(ymin=value-se,ymax=value+se),width=0.1)
	p <- p +geom_line()+geom_point() + ggtitle(curgene)
	p <- p +xlab("Post-conception week (BrainSpan)")+ylab("Expression")
	p <- p + ylim(c(3,11))
	# add dev milestone references
	p <- BS_addRefs(p)
	p <- p + coord_trans(x="log10")
	plotList[[curgene]] <- p
}
	pdf(sprintf("%s/BrainSpan_indivGenes_%s.pdf",outDir,dt),width=11,height=9)
	tryCatch({
		nr=3;nc=1;
	  for (k in seq(1,length(plotList),nr*nc)) {
	    print(multiplot(plotlist=plotList[k:(k+((nr*nc)-1))],
	        layout=matrix(1:(nr*nc),nrow=nr,ncol=nc)))
	  }
	},error=function(ex){print(ex)},finally={dev.off()})
}
}

# -------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------

#' convert text age into days from conception
BS_convertAge <- function(age) {
	age <- as.character(age)
	age <- sub("age: ","",age)
	new_age <- rep(NA,length(age))
	i1 <- grep("PCW",age)			# sel age
	new_age[i1] <- as.integer(sub(" PCW","",age[i1]))*7		# week->days
	i2 <- grep("PMonth",age)	# first year of life
	new_age[i2] <- as.integer(sub(" PMonth","",age[i2]))*30 # month->days
	new_age[i2] <- new_age[i2] + (40*7) # add days for 40w gestation
	i3 <- grep("Y",age)				# rest of cycle
	new_age[i3] <- as.integer(sub(" Y","",age[i3]))*365 # year->days
	new_age[i3] <- new_age[i3] + (40*7) # add days for gestation
	
	return(new_age)
	
}

#' add BSpan dev milestone labels + reference lines
#' p (ggplot2 object)
BS_addRefs <- function(p,yval=105) {
# add reference lines and life markers
first_trim	<- 12*7
second_trim <- 24*7
third_trim	<- 40*7
first_year	<- (40+52)*7
adol				<- third_trim + (12*365)
adult				<- third_trim + (21*365)
#p <- p+geom_vline(xintercept=c(first_year+ (1:12)*365),lty=2)
p <- p+geom_vline(xintercept=c(first_trim,second_trim,third_trim,
															first_year,adol,adult),lty=2,colour="black")
p <- p+annotate("text",x=first_trim,y=yval,label="1st trim.",hjust=1,fontface=4,colour="blue")
p <- p+annotate("text",x=second_trim,y=yval,label="2nd trim.",hjust=1,fontface=4,colour="blue")
p <- p+annotate("text",x=third_trim,y=yval,label="3rd trim.",hjust=1,fontface=4,colour="blue")
p <- p+annotate("text",x=first_year,y=yval,label="Pnatal 1yr",hjust=1,fontface=4,colour="blue")
p <- p+annotate("text",x=adol,y=yval,label="Child",hjust=1,fontface=4,colour="blue")
p <- p+annotate("text",x=adult,y=yval,label="Teen",hjust=1,fontface=4,colour="blue")
#p <- p+scale_x_discrete(breaks=c(first_trim,second_trim,third_trim,first_year,adol,adult))
p
}

