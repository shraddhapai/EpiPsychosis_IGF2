#' plot betas by group, for selected probes
#' 
#' @param mset (GenomicRatioSet)
#' @param group (char) which group variable to separate probes
#' @param baseline_group (char) which class in 'group' is the baseline?
#' function draws a reference line at the median for that value
#' @param group2 (char; optional) split boxplots within a group using a 
#' second conditional variable. e.g. group="DX" and group2="SEX",
#' will show a primary split by diagnosis and secondary by sex. If NULL
#' not used
#' @param selProbes (char) probes to plot
#' @param rmLegend (logical) if TRUE does not show any legends. May be
#' used when plotting several in one page for a manuscript.
#' @param plot_nr (integer) number rows in plot
#' @param plot_nc (integer) number columns in plot
#' @export
#' @import ggplot2
dmp_plotProbes_groups <- function(mset,group=NULL,baseline_group=NULL,
					  group2=NULL,selProbes,rmLegend=FALSE,plot_nr=3,plot_nc=4) {
	mset <- mset[which(rownames(mset)%in% selProbes),]
	betas <- getBeta(mset)
	pd <- pData(mset)[,group]

	if (!is.null(group2)){
		pd2 <- pData(mset)[,group2]
	}

	if (!is.factor(pd)) pd <- factor(pd)

	loc <- getLocations(mset)
	#par(mfrow=c(4,3),bty='n',las=1)
	pset <- list()
	for (i in 1:nrow(betas)) {
		g2str <- ""
		if (is.null(group2)) {
			df <- data.frame(GROUP=pd, betas=betas[i,])
			p <- ggplot(df, aes(GROUP,betas)) 
			p <- p + geom_boxplot()
		} else {
			df <- data.frame(GROUP=pd, betas=betas[i,],GROUP2=pd2)
			p <- ggplot(df, aes(GROUP,betas)) 
			p <- p + geom_boxplot(aes(colour=GROUP2))
			g2str <- group2
		}
		ttl <- sprintf("%s\n%s;%s",names(loc)[i],group,g2str)

		md <- median(df$betas[which(df$GROUP%in%baseline_group)])
		p <- p + geom_hline(yintercept=md)

		p <- p + ggtitle(ttl) # + ylim(0,1)
		print(ttl)

		if (rmLegend) {
			p <- p+ theme(legend.position="none")
		}
																			 		pset[[i]] <- p
	}

	nr=plot_nr;nc=plot_nc;
	for (k in seq(1,length(pset),nr*nc)) {
		print(multiplot(plotlist=pset[k:(k+((nr*nc)-1))], 
				layout=matrix(1:(nr*nc),nrow=nr,ncol=nc)))
	}
}
