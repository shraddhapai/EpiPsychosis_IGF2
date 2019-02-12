# pathway analysis mouse results

gmtFile <- "/home/shraddhapai/Epigenetics/NARSAD/anno/Mouse_GOBP_AllPathways_no_GO_iea_October_01_2017_symbol.gmt"
gsea <- "/home/shraddhapai/software/gsea-3.0.jar"

sigFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/MouseRnaseq2/diffEx_180821/gpSTR/diffEx.txt"

dat <- read.delim(sigFile,sep="\t",h=T,as.is=T)
dat$score <- -log10(dat$PValue)*sign(dat$logFC)
dat <- dat[,c("gene.name","score")]
dat <- dat[order(dat[,2]),]
rnkFile <- sprintf("%s/diffEx.Pvalues.GSEA_input.rnk",dirname(sigFile))
write.table(dat,file=rnkFile,sep="\t",col=T,row=F,quote=F)
outDir <- dirname(sigFile)
gsea_cmd <- paste("java -cp ", gsea, " -Xmx8G xtools.gsea.GseaPreranked -gmx ", gmtFile, " -norm meandiv -nperm 1000 -rnk ", rnkFile, " -scoring_scheme classic -rpt_label mouseSTR_RNAseq -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 200 -set_min 10 -zip_report false -out ", outDir, " -gui false", sep="")
system(gsea_cmd)

