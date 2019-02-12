#' compare assigned and predicted sex
#'
#' @details Column should have values "Female" or "Male". Should also have
#' a "Sample_Name" column
#' @param mset (MethylSet)
#' @param pheno (char) column name. 
#' @return (char) vector of Sample_Names where assigned and predicted sex
#' don't match
checkSex <- function(mset,pheno="SEX") {

	pd <- pData(mset)
	sex <- pd[,pheno]; 
	sex[which(pd[,pheno] %in% "Female")] <- "F"; 
	sex[which(pd[,pheno] %in% "Male")] <- "M"

	print(table(pd$predictedSex,sex))
	idx <- which(sex!=pd$predictedSex & sex %in% c("F","M"))
	if (length(idx)>0) { 
		cat(sprintf("\tSamples with mismatched sex\n"))
		warnSamp <- pd$Sample_Name[idx]
	} else {
		cat("Woo hoo! No samples with mismatched sex!\n")
		warnSamp <- NULL
	}
		
	return(warnSamp)
}
