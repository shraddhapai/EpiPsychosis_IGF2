#' explore out of bound control probes 
#'
#' @details this script gets the address of the offending probes for 
#' communication with Illumina
#' @param rgSet (RGSet) raw read-in values
#' @param xlim (numeric, length 2) acceptable range of intensities
#' @param sampNames (char) sample names for reference
#' @import minfi
#' @export
controlProbes_oob <- function(rgSet,xlim=c(5,17),sampNames=NULL) {
    r <- getRed(rgSet)
    g <- getGreen(rgSet)
		
control <- c("BISULFITE CONVERSION I","BISULFITE CONVERSION II","SPECIFICITY I")
for (controlType in control) {
	ctrlAddress <- getControlAddress(rgSet, controlType = controlType)
	ctlWide <- log2(r[ctrlAddress, ])
	if (!is.null(sampNames)) 
		colnames(ctlWide) <- sampNames
	ctlR <- melt(ctlWide, varnames = c("address", "sample"))
	ctlWide <- log2(g[ctrlAddress, ])
	if (!is.null(sampNames)) 
	colnames(ctlWide) <- sampNames
	ctlG <- melt(ctlWide, varnames = c("address", "sample"))
	ctl <- rbind(cbind(channel = "Red", ctlR), cbind(channel = "Green", 
	ctlG))
	
	if (any((ctl$value < xlim[1]) | (ctl$value > xlim[2])))  {
		message("Warning: ", controlType, " probes outside plot range")
		ctl <- ctl[order(ctl$address),]
		idx <- which(ctl$value < xlim[1] | ctl$value>xlim[2])
		print(ctl[idx,])
		cat("-------------\n")
	} else {
		cat("\tNo problems\n")
	}

}
}
