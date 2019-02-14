require(minfi)
load('~/Desktop/paper writing/SCZ/comb-p/CaseControlEPIC_CLEAN_161221.Rdata')

locs <- getLocations(MSet.genome)
loc <- data.frame(names(locs), data.frame(locs))
colnames(loc) <- c("name", "seqnames", "start", "end", "width", "strand")


dat_agesexPMI_CEU_PC12_noSlide_171129 <- read.delim("agesexPMI_CEU_PC12_noSlide_171129.txt")
dat_agesexPMI_CEU_PC12_noSlide_171129 <- data.frame(rownames(dat_agesexPMI_CEU_PC12_noSlide_171129), dat_agesexPMI_CEU_PC12_noSlide_171129)
colnames(dat_agesexPMI_CEU_PC12_noSlide_171129) <- c("name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
head(dat_agesexPMI_CEU_PC12_noSlide_171129)

library(data.table)
keys <- c("name")

dat_agesexPMI_CEU_PC12_noSlide_171129_table <- data.table(dat_agesexPMI_CEU_PC12_noSlide_171129, key=keys)
loc_table <-  data.table(loc, key=keys)
cgInfor_agesexPMI_CEU_PC12_noSlide_171129 <- loc_table[dat_agesexPMI_CEU_PC12_noSlide_171129_table, nomatch=0]
cgInfor_agesexPMI_CEU_PC12_noSlide_171129 <- data.frame(cgInfor_agesexPMI_CEU_PC12_noSlide_171129)
length(dat_agesexPMI_CEU_PC12_noSlide_171129[,1])
length(cgInfor_agesexPMI_CEU_PC12_noSlide_171129[,1])

head(cgInfor_agesexPMI_CEU_PC12_noSlide_171129)
write.table(cgInfor_agesexPMI_CEU_PC12_noSlide_171129, "cgInfor_agesexPMI_CEU_PC12_noSlide_171129.txt", row.names=FALSE, quote=FALSE)

cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final <- data.frame(cgInfor_agesexPMI_CEU_PC12_noSlide_171129$seqnames, cgInfor_agesexPMI_CEU_PC12_noSlide_171129$start, cgInfor_agesexPMI_CEU_PC12_noSlide_171129$end, cgInfor_agesexPMI_CEU_PC12_noSlide_171129$P.Value, cgInfor_agesexPMI_CEU_PC12_noSlide_171129$t)
head(cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final)
write.table(cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final, "cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final.txt", row.names=FALSE, quote=FALSE)

write.table(cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final, "cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final.bed", row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")


sort -k1,1 -k2,2n -o /Users/peipei.li/Desktop/paper\ writing/SCZ/comb-p/171129/cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final_sort.bed /Users/peipei.li/Desktop/paper\ writing/SCZ/comb-p/171129/cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final.bed

comb-p pipeline -c 4 --seed 0.01 --dist 500 --region-filter-p 0.01 -p 171129/agesexPMI_CEU_PC12_noSlide_171129/out_agesexPMI_CEU_PC12_noSlide_171129 --annotate hg19 /Users/peipei.li/Desktop/paper\ writing/SCZ/comb-p/171129/cgInfor_agesexPMI_CEU_PC12_noSlide_171129_final_sort.bed

