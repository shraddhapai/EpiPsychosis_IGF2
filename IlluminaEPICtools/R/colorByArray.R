#' color phenotype by position on array
#'
#' @param v (factor) variable of interest, one value per array
#' @param arrayNum (integer) which array will sample be on? same length as v
#' @param palName (char) RColorBrewer palette for colour-coding
#' @param ttl (char) plot title
#' @param nr (integer) num rows per array
#' @param nc (integer) num cols per array
#' @param arrayPos (integer) location within the array. if NULL, uses 
#' order in v
#' @import grid
#' @export
colorByArray <- function(v, arrayNum,palName="Dark2",ttl="Array layout",
		nr=4,nc=2,arrayPos=NULL) {
	pal <- suppressMessages(brewer.pal(name=palName,n=length(levels(v))))
	maxArray <- max(arrayNum)

	cdim <- 0.05		# width of cell
	array_gap <- 0.1	# space between arrays minus cdim
	spos <- 0			# start x pos
	ypos <- 0			# start y pos
	ldim <- 0.03		# legend cell dim

	# opening the graphic device and
	# setting up a viewport with borders:
	grid.newpage()
	vp1 <- viewport(x=0.1,y=0.1,width=0.8,height=0.8,name="vp1",
					just=c("left","bottom"))
	grid.text(ttl, x=0.5,y=1,vp=vp1)
	# legend
	lv <- levels(v); 
	n <- length(lv)
	legd <- data.frame(x=rep(0.8,n),y=seq(0.94,0.94-((n-1)*ldim),-ldim))
	grid.rect(x=legd$x,y=legd$y,gp=gpar(col=1,fill=pal[1:n]),
			  width=ldim,height=ldim, vp=vp1)
	grid.text(lv,x=0.8+(2.3*ldim),y=legd$y,vp=vp1,hjust="right")

	for (k in 1:maxArray) {
		cat(sprintf("Array %i: (%1.2f, %1.2f)\n", k, spos,ypos))
		subv <- v[which(arrayNum %in% k)]
		if (!is.null(arrayPos)) {
			subpos	<- arrayPos[which(arrayNum %in% k)]
		print(subpos)
			subv	<- subv[subpos]; 
		}
		if (length(subv) != nr*nc) {
			cat(sprintf("\tNum samples (%i) does not equal num slots (%i)\n",
				length(subv),nr*nc))
		}
		dat <- data.frame(x = rep(seq(spos,spos+(cdim*(nc-1)),cdim), nr), 
					  y = rep(seq(ypos,ypos+(cdim*(nr-1)),cdim), each = nc),
					  col=pal[as.integer(subv)])

		# plotting rectangles using x/y positions
		grid.rect(x=dat$x,y=dat$y,
				  height=cdim,width=cdim,
				  hjust=0,vjust=0,vp=vp1,
		          gp=gpar(col=1, fill=as.character(dat$col)))
		# array name
		grid.text(sprintf("array %i",k), x=spos,y=ypos+((nr+0.5)*cdim),
				  vp=vp1,
				  just=c("left","bottom"))
				  
		spos <- spos + (cdim*(nc-1)) + array_gap
		if (spos > 0.9) { # move up a row
			spos <- 0; 
			ypos <- ypos + (cdim*nr) + array_gap
		}
	}
}
