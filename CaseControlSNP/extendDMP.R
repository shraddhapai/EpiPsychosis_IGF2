#' extend dmp regions for cis-eqtl analysis
args <- commandArgs(TRUE)
dmpFile <- args[1]

require(GenomicRanges)
dat <- read.delim(dmpFile,sep="\t",h=T,as.is=T)
dat <- subset(dat,z_sidak_p<0.05)
gr <- GRanges(dat[,1],IRanges(dat[,2],dat[,3]))
# get 1mb window centered on dmp
gr <- resize(gr, width=1000000,fix="center")
df <- as.data.frame(gr)
options(scipen=10)
write.table(df[,c(1,2,3)],
    file=sprintf("%s/dmpQ0.05.1Mb.win.txt",dirname(dmpFile)),
	sep="\t",col=F,row=F,quote=F)

