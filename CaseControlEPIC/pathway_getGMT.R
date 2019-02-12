#' quick script to read pathway ORA file, fix the p/Q values (weren't
#' correctly computed due to bug in dmp_pathwayORA) and generate GMT file

require(IlluminaEPICtools)
require(minfi)

rootDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398"
inFile <- sprintf("%s/preprocessing/quantile/CaseControlEPIC_MSetgenome_Qtl_160923.Rdata",rootDir)
dmpFile <- sprintf("%s/dmp_DX_QTL/dmp.blockTechReps.dropSNPs.dmp_factor(DX)control.topTable.160929.txt",rootDir)

pResFile <- sprintf("%s/dmp_DX_QTL/dmp_pathways_stats_pvalue.0.01.161014.txt", rootDir)
annoFile <- "Annotations.txt"

# -----------------
cat("* load MSet.genome")
print(system.time(load(inFile)))

cat("* load dmp results and subset foreground probes\n")
dmp <- read.delim(dmpFile,sep="\t",h=T,as.is=T)
dmp <- subset(dmp, P.Value < 0.005)
fg_GR <- getLocations(MSet.genome[rownames(dmp)])

cat("* compile genes and pathways")
annoFiles <- read.delim(annoFile,sep="\t",h=T,as.is=T)
genes <- read.delim(annoFiles[which(annoFiles[,1]=="Genes_Symbols"),3],
			sep="\t",h=F,as.is=T)
gene_GR <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),name=genes[,4],
	strand=genes[,6])
pathwayFile <- annoFiles[which(annoFiles[,1]=="Pathways"),3]

cat("* Get pathway results\n")
pRes <- read.delim(pResFile,sep="\t",h=T,as.is=T)
pRes$Hypergeom_Q <- p.adjust(pRes$Hypergeom_p,method="BH")
pRes$Binom_Q <- p.adjust(pRes$Binom_p,method="BH")

FDRcutoff <- 0.1
selPathways <- pRes$Pathway[which(pRes$Hypergeom_Q < FDRcutoff)]
outFile <- sprintf("%s.pathwayp0.005.Q%1.2f.gmt",sub(".txt","",dmpFile), 
				   FDRcutoff)
dmp_pathwayORA_writeGMTFromRes(pathwayFile,selPathways,fg_GR,
	gene_GR, geneDomain,outFile)


