#' extend dmp hits 1mb

dmpFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/age_sex_c1c2_slide/dmp.blockTechReps.dropSNPs.dmp_factor(DX)control.topTable.161205.FDR0.30.locations.txt"

dmp <- read.delim(dmpFile,sep="\t",h=T,as.is=T)
dmp <- subset(dmp,qval < 0.05)
require(GenomicRanges)

gr <- GRanges(dmp[,2],IRanges(dmp[,3],dmp[,3]),name=dmp[,1])
gr <- resize(gr, width=1000000,fix="center")
df <- as.data.frame(gr)
options(scipen=10)
outFile <- sprintf("%s.1Mbwins.txt",dmpFile)
write.table(df[,c(1:3,6)],file=outFile,sep="\t",col=T,row=F,quote=F)


