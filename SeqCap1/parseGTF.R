#' parse GTF file to return gene definition in table
#' @param gFile (char) path to gtf file
#' @param header (logical) has header row
#' @return (data.frame) one row per record. chrom. start. end. strand. ENGSID.
#' gene_name
parseGTF <- function(gFile,header=F) {

fName <- file(gFile,"r")
out <- list()
ctr <-1
repeat {
    currLine <- scan(gFile,what="character",nlines=1,quiet=TRUE,sep="\t")
    if (length(currLine)==0) break
	info <- unlist(strsplit(currLine,";"))
	id <- sub("gene_id ","",info[grep("gene_id",info)])
	id <- gsub("\"", "",id)
	name <- sub("gene_name ","",info[grep("gene_name",info)])
	name <- gsub("\"", "",name)
	out[[ctr]]<- c(currLine[c(1,4,5,7)],id,name)
	ctr <- ctr+1
}
out <- do.call("rbind",out)
colnames(out) <- c("chr","start","end","strand","gene_id","gene_name")
out <- as.data.frame(out)
out

}
