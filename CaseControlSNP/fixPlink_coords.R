
# plink coords for SNP data are incorrect in the SUS19399.bim file
# obtained from TCAG.
# spot checks of rs4601751 and rs12575480 confirm this.
# This script updates the coordinates for the SNPs using the manifest
# file from the illumina support website. 

# The correct SNP coordinates, as verified against the ucsc browser,
# are found in the bpmap file available on Illumina's website.
# http://support.illumina.com/array/array_kits/infinium-psycharray-beadchip-kit.html
# http://support.illumina.com/downloads/infinium-psycharray-v1-1-product-files.html
# File: ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/infinium-psycharray/v1-1/infinium-psycharray-24-v1-1-a1-csv.zip
# Last updated 20 Nov 2016

# ------------------------------------------
# SNP coords were extracted from that file with the following code
# cat InfiniumPsychArray-24v1-1_A2.csv | cut -f1-4,7-11 -d',' > tmp.txt
# cat tmp.txt | awk -F"," '{print "chr"$8"\t"$9"\t"$9"\t"$2}' > PsychArray_24.bed
# these coordinates are known to be hg38 coordinates.

bedFile <- "/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399_New/hg38/PsychArray_24.bed"
bimFile <-  "/scratch/g/gbader/spai/NARSAD2014/output_files/SUS19399_New/hg38/SUS19399.bim"

cat("reading bim\n")
bim <- read.table(bimFile,h=F,as.is=T)
cat("Reading bed\n")
coords <- read.delim(bedFile,sep="\t",h=T,as.is=T)

colnames(coords) <- c("chrom","start","end","name")
colnames(bim) <- c("chrom","name","gencoord","loccoord","a1","a2")

midx <- match(bim$name,coords$name)
# TRUE - sp manual check
all.equal(coords$name[midx], bim$name)
coords <- coords[midx,]

bim[,4] <- coords$start
options(scipen=10)
outF <- sprintf("%s.fixed",bimFile)
write.table(bim,file=outF,sep="\t",col=F,row=F,quote=F)





