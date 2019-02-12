#' extract SNPs that overlap the vicinity of case/control DMPs and nearby
#' regulatory regions
rm(list=ls())

dmpRegionFile <- "/scratch/g/gbader/spai/NARSAD2014/output_files/CaseControlEPIC/dmp_DX_QTL/dmpCaseControl_Genes_161206.txt"
outDir <- "/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399/meQTL"

dat <- read.delim(dmpRegionFile,h=F,as.is=T)
require(GenomicRanges)
gr <- GRanges(dat[,1],IRanges(dat[,2],dat[,3]),name=dat[,4])
gr <- resize(gr,width=width(gr)+18000,fix="center")
df <- as.data.frame(gr)
head(df)
options(scipen=10)
outFile <- sprintf("%s/%s.20K.txt",outDir,
        sub(".txt","",basename(dmpRegionFile)))
print(outFile)
write.table(df[,c(1,2,3,6)],
    file=outFile,
    sep="\t",col=F,row=F,quote=F)
