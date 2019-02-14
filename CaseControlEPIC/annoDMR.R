#' convert dmr results to bed file for comb-p
#' compute bacon-adjusted pvalues 

dmrFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/agesexPMI_PC12_171129/dmp..dmp_DXcase.topTable.171129.txt"
datFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata"

require(minfi)
require(bacon)
load(datFile)
locs <- as.data.frame(getLocations(MSet.genome))
locs$probeID <- rownames(locs)

dmr <- read.delim(dmrFile,sep="\t",h=T,as.is=T)
bc <- bacon(dmr$t)
cat(sprintf("Inflation rate = %1.2f\n", inflation(bc)))
dmr$corrP <- pval(bc)
dmr$probeID <- rownames(dmr)

x <- merge(x=locs,y=dmr, by="probeID") 
x <- x[,c("seqnames","start","end","corrP")]

write.table(x,file=sprintf("%s_corrP.bed",dmrFile),sep="\t",col=F,row=F,quote=F) 

