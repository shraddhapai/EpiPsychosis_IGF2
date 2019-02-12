require(Rsamtools)
require(doParallel)

#' returns methylation for a set of GRanges from an input tabix file
#' 
#' @details The tabix file must be in the output format of methratio.py from the 
#' bsmap package. 
#' chr	pos	strand	context	ratio	eff_CT_count	C_count	CT_count	
#' rev_G_count	rev_GA_count	CI_lower	CI_upper
#' @param tbxFile (character) tabix file with columns 
#' @param gr (GenomicRanges) a set of ranges to retrieve methylation for
#' @param context (character) currently for single context. Looks at the third
#' and fourth letter of the context (e.g. GGCGA would be "CG" but "GACCT" would
#' be CC"
#' @param unit (character) one of [pctM | MU | Cov]
#' @param STR_POS, CONTEXT_POS, U_POS, COV_POS (integer) 
#' column of tabix file with strand, context
#' U position, (M+U) position. 
#' @param numFields (integer) num fields in a tabix record.
#' @param splitStrands (logical) If TRUE, does not pool reads from W and C 
#' strands
#' 	if FALSE, does pool
#' @param getBaseLevel (logical): does not aggregate for all bases of 
#' a locus. 
#' 	When TRUE, ignores splitStrands, and returns reads as retrieved from
#' 	tabix file
#' @param printRecords (logical) prints all records for each input GRanges to 
#' standard out. Intended for debugging purposes only.
#' @param base_minCvg (integer) for locus-level aggregation, apply coverage
#' filter for bases to be included in aggregation. If splitStrands=TRUE,
#' filter is applied for each strand. Set to NULL for no filter.
#' @param numCores (integer) number of cores for parallel processing
#' balance with known memory requirements for each record to process
getRec_GRanges <- function(tbxFile, gr, context="CG",unit="pctM",
	STR_POS=3,CONTEXT_POS=4,U_POS=7,MU_POS=8,numFields=12,splitStrands=TRUE,
    printRecords=FALSE, base_minCvg=1L,numCores=1L,verbose=FALSE
) {

hdr <-c("chr","pos","strand","context","ratio","eff_CT_count","C_count","CT_count","rev_G_count","rev_GA_count","CI_lower","CI_upper")

if (any(strand(gr)%in% c("+","-"))) {
		cat("gr should have a wildcard strand. Forcing.\n")
		strand(gr) <- "*"
}

f		<- Rsamtools::TabixFile(tbxFile)
snames	<- headerTabix(f)$seqnames
sno_idx	<- which(!seqnames(gr)%in% snames)
gnames	<- gr$name

if (any(sno_idx)) {
		cat("ranges not in seqnames(tabix) detected\n")
		gr	<- gr[-sno_idx]
}
    if (verbose) cat("\tFetching records\n")
    t0 <- Sys.time()
res	<- scanTabix(f, param=gr)
    cat(sprintf("\t%i records (%1.3f min)\n", length(res),
			(Sys.time()-t0)/60))

    #cat("\tGetting unit-wise measures\n")
    t0 <- Sys.time()
out <- mclapply(res, function(m) {
	emptyFlag	<- FALSE
	if (length(m)==0) emptyFlag	<- TRUE
	else {
		rec <- unlist(strsplit(m,"\t"))
	tryCatch({
		rec	<- matrix(rec,byrow=TRUE,ncol=numFields)
	},error=function(ex){
		print(ex)
		browser()
	})

	if (printRecords) print(rec)
	rec[,CONTEXT_POS] <- substr(rec[,CONTEXT_POS],3,4)
	valid_ctxt <- NA
	if (context == "CG") valid_ctxt <- c("CG","GC")
	else { stop("context other than CG not implemented\n")}
		idx <- which(rec[,CONTEXT_POS]%in% valid_ctxt)
		if (verbose) {
			cat(sprintf("\t%i of %i are %s context\n", length(idx),
				nrow(rec),context))
		}
	if (length(idx)<1) emptyFlag <- TRUE
	else {
		rec <- data.frame(rec[idx,,drop=FALSE])
		rec[,U_POS] <- as.integer(as.character(rec[,U_POS]))
		rec[,MU_POS] <- as.integer(as.character(rec[,MU_POS]))
		colnames(rec) <- hdr
	}
}
if (emptyFlag) val <- NA else val <- rec
	return(val)
},mc.cores=numCores)

# create placeholders for input ranges that were from
# seqnames that were not present in the tabix file
if (any(sno_idx)) {
	cat("this block of code hasn't been implemented and needs to be checked"); browser()
	n1		<- length(gr); n2 <- length(sno_idx)
	pres	<- setdiff(1:(n1+n2),sno_idx)
	out2	<- matrix(NA,nrow=n1+n2,ncol=ncol(out))
	out2[pres,] 			<- out
	rownames(out2)[pres] 	<- names(out2)
	rownames(out2)[sno_idx]	<- gnames[sno_idx]
	browser()
}
out
}

