#-----------------------------------------------------------------------------------------
#-------------------------------------- SCZ RNAseq ---------------------------------------
#-----------------------------------------------------------------------------------------

#--------- 2017.11.27
#--------- Untar Ungzip
qsubs 1 24 'tar -zxvf LABV_20170302_RNA.tar.gz' untar_ungzip


#--------- 2017.11.28
#--------- Gencode
#- Human Gencode_GRCh38_p10
https://www.gencodegenes.org/releases/27.html
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v27.primary_assembly.annotation.gtf.gz
#- Mouse Gencode_GRCm38_p5
https://www.gencodegenes.org/mouse_releases/current.html
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.primary_assembly.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
gunzip gencode.vM15.primary_assembly.annotation.gtf.gz

#--------- STAR INDEX
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/star/2.5.3a/bin/Linux_x86_64
STAR --help
mkdir star_index_75bp

qsubs 10 24 '\
STAR --runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir star_index_75bp \
	--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
	--sjdbGTFfile gencode.v27.primary_assembly.annotation.gtf \
	--sjdbOverhang 75 \
' star_index_75bp

qsubs 10 24 '\
STAR --runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir star_index_75bp \
	--genomeFastaFiles GRCm38.primary_assembly.genome.fa \
	--sjdbGTFfile gencode.vM15.primary_assembly.annotation.gtf \
	--sjdbOverhang 75 \
' star_index_75bp

#--------- RSEM INDEX
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/rsem/RSEM/bin
rsem-prepare-reference --help
mkdir rsem_index

qsubs 6 24 '\
rsem-prepare-reference --num-threads 6 \
	--gtf gencode.v27.primary_assembly.annotation.gtf  \
	GRCh38.primary_assembly.genome.fa \
	rsem_index/rsem \
' rsem_index

qsubs 6 24 '\
rsem-prepare-reference --num-threads 6 \
	--gtf gencode.vM15.primary_assembly.annotation.gtf  \
	GRCm38.primary_assembly.genome.fa \
	rsem_index/rsem \
' rsem_index

#--------- FASTQ MERGE
mkdir fastq
for i in `cat names` ; do 
qsubs 1 24 '\
cat run[1-4]/'${i}'_L000_R1_001.fastq.gz > fastq/'${i}'_R1.fastq.gz
cat run[1-4]/'${i}'_L000_R2_001.fastq.gz > fastq/'${i}'_R2.fastq.gz \
' fastq_merge
done

#while read line ; do 
#	mv $(awk '{print $1"_R2.fastq.gz", $2"_R2.fastq.gz"}' <<<$line) 
#done < ../names

#--------- FASTQC
mkdir fastqc
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/fastqc/0.11.5/
fastqc --help

for i in `cat names` ; do 
qsubs 1 24 '\
fastqc --outdir fastqc \
	--format fastq \
	fastq/'${i}'_R1.fastq.gz \
	fastq/'${i}'_R2.fastq.gz \
' fastqc_pretrim
done

#--------- TRIMGALORE
mkdir trimgalore
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/fastqc/0.11.5/
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/trimgalore/0.4.4/
trim_galore --help

for i in `cat names` ; do 
qsubs 1 24 '\
trim_galore --quality 20 \
	--phred33 \
	--fastqc \
	--length 50 \
	--output_dir trimgalore \
	--paired \
	fastq/'${i}'_R1.fastq.gz \
	fastq/'${i}'_R2.fastq.gz \
' trimgalore
done


#--------- STAR ALIGN
mkdir star
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/star/2.5.3a/bin/Linux_x86_64
STAR --help

for i in `cat names | sed -n '1,34p'` ; do 
qsubs 4 24 '\
STAR --runMode alignReads \
	--runThreadN 4 \
	--genomeDir ~/labrie-primary/Bioinformatics_core/genomes/Human/Gencode_GRCh38_p10/star_index_75bp \
	--genomeLoad NoSharedMemory \
	--readFilesIn trimgalore/'${i}'_R1_val_1.fq.gz trimgalore/'${i}'_R2_val_2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix star/'${i}'_ \
	--outReadsUnmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outFilterMismatchNoverLmax 0.1 \
	--quantMode TranscriptomeSAM \
' star_aligning
done

for i in `cat names | sed -n '35,39p'` ; do 
qsubs 4 24 '\
STAR --runMode alignReads \
	--runThreadN 4 \
	--genomeDir ~/labrie-primary/Bioinformatics_core/genomes/Mouse/Gencode_GRCm38_p5/star_index_75bp \
	--genomeLoad NoSharedMemory \
	--readFilesIn trimgalore/'${i}'_R1_val_1.fq.gz trimgalore/'${i}'_R2_val_2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix star/'${i}'_ \
	--outReadsUnmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outFilterMismatchNoverLmax 0.1 \
	--quantMode TranscriptomeSAM \
' star_aligning2
done



#--------- RSEM Counting
mkdir rsem
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/rsem/RSEM/bin
rsem-calculate-expression --help

for i in `cat names | sed -n '1,34p'` ; do 
qsubs 4 24 '\
rsem-calculate-expression \
	--paired-end \
	--num-threads 4 \
	--bam \
	--no-bam-output \
	star/'${i}'_Aligned.toTranscriptome.out.bam \
	~/labrie-primary/Bioinformatics_core/genomes/Human/Gencode_GRCh38_p10/rsem_index/rsem \
	rsem/'${i}' \
' rsem_transcriptome
done

for i in `cat names | sed -n '35,39p'` ; do 
qsubs 4 24 '\
rsem-calculate-expression \
	--paired-end \
	--num-threads 4 \
	--bam \
	--no-bam-output \
	star/'${i}'_Aligned.toTranscriptome.out.bam \
	~/labrie-primary/Bioinformatics_core/genomes/Mouse/Gencode_GRCm38_p5/rsem_index/rsem \
	rsem/'${i}' \
' rsem_transcriptome
done

#--------- MULTIQC
python -V 
qsubs 4 24 'multiqc .' multiqc_reports












#-----------------------------------------------------------------------------------------
#--------------------------- SCZ RNAseq RSEM EdgeR Cibersort -----------------------------
#-----------------------------------------------------------------------------------------
qrsh4
mkdir expected_count && cd expected_count
for i in `cat ../../names` ; do
awk -v OFS='\t' '{print $1, $5}' ../../rsem/${i}.genes.results | sed -e 's/\..*\t/\t/' > ${i}
done

mkdir rpkm_count && cd rpkm_count
for i in `cat ../../names` ; do
awk -v OFS='\t' '{print $1, $7}' ../../rsem/${i}.genes.results | sed -e 's/\..*\t/\t/' > ${i}
done

for i in `cat ../../names` ; do 
qsubr 1 24 '
library(readr)
sample <- read_delim("'${i}'","\t",escape_double=FALSE,trim_ws=TRUE)
sample <- sample[!duplicated(sample),]
write.table(sample,file="'${i}'_edit",sep="\t",col=T,row=F,quote=F)
' r_duplicate
done

#--------- R EdgeR
library(readr)
covar <- read_delim("covar.tsv","\t",escape_double=FALSE,trim_ws=TRUE)

#- Merge Count data
library(edgeR)
counts <- paste(covar$ID,"_edit",sep="")
counts <- readDGE(counts,path="expected_count",columns=c(1,2),labels=covar$ID,header=TRUE)
counts <- counts$counts
write.csv(counts,file="SCZ_RNAseq_STAR_RSEM_counts.csv")

rpkm <- paste(covar$ID,"_edit",sep="")
rpkm <- readDGE(rpkm,path="rpkm_count",columns=c(1,2),labels=covar$ID,header=TRUE)
rpkm <- rpkm$counts

#- DGEList
library(edgeR)
y <- DGEList(counts, group = covar$DISEASE)

#- Filtering
keep <- rowSums(cpm(y)>1) >= 34
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)		#[1] 15322    34

#- Normalization
y <- calcNormFactors(y,method="TMM")
y$samples

#- CPM
cpm <- cpm(y,normalized.lib.sizes=TRUE, log=FALSE)
write.csv(cpm,file="cpm.csv")
logcpm <- cpm(y,normalized.lib.sizes=TRUE, log=TRUE)
write.csv(logcpm, file="logcpm.csv")

#--------- Dendogram
library(dendextend)
dend <- as.dendrogram(hclust(dist(t(cpm))))
par(cex=1.2, mar=c(4, 4, 4, 0))
labels_cex(dend) <- 0.4
plot(dend, type = "rectangle", ylab = "Height", main = "SCZ RNAseq STAR RSEM CPM")
#-Colors Labels
colors <- c("green3","orange3","red3")[factor(covar$DISEASE)]
labels_colors(dend) <- colors[order.dendrogram(dend)]
plot(dend, type = "rectangle", ylab = "Height", main = "SCZ RNAseq STAR RSEM CPM")
#-Colour Bars
library(RColorBrewer)
display.brewer.all()
the_bars <- cbind(
	heat.colors(39)[covar$PMI],
	terrain.colors(77)[covar$AGE],
	topo.colors(19)[covar$RIN],
	brewer.pal(4,"Dark2")[as.factor(covar$`Extraction Date`)],
	c("green3","orange3","red3")[factor(covar$DISEASE)])
par(mar = c(7,5,5,1))
labels_cex(dend) <- 0.8
plot(dend, type = "rectangle", ylab = "Height", main = "SCZ RNAseq STAR RSEM CPM")
colored_bars(colors=the_bars,dend=dend,rowLabels=c("PMI","AGE","RIN","DATE","DISEASE"),
	cex.rowLabels=0.7,y_shift=-16000,text_shift=0)
#--------- MDSplot
colors <- c("green3","orange3","red3")[factor(covar$DISEASE)]
plotMDS(y,col=colors,main="SCZ RNAseq STAR RSEM CPM")
legend("topright",fill=c("green3","orange3","red3"),legend=c("CTRL","BIPOL","SCZ"),cex=1)
#--------- PCAplot
pca <- prcomp(t(cpm), scale = TRUE)
library("ggbiplot")
ggbiplot(pca,choices=1:2,obs.scale=1,var.scale=1,groups=covar$DISEASE,ellipse=TRUE,
	circle=TRUE,var.axes=F,labels=covar$DISEASE) +
	scale_color_discrete(name="SCZ RNAseq STAR RSEM CPM") +
	theme(legend.direction = 'horizontal', legend.position = 'top')
#--------- Correlation
library(corrplot)
corr <- cor(cpm,use="complete.obs", method ="pearson")
col1 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
                           "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(corr, method="circle", col = col1(200), title = "SCZ RNAseq STAR RSEM CPM",
    mar = c(1, 1, 2, 1), order="hclust", hclust.method = c("complete"), addrect=5, tl.cex = 0.7,
    tl.col = "black")






#--------- CIBERSORT
#- Sort cibersort signature from Qianhui_Yu paper
library(readxl)
cibersort_signatures <- read_excel("cibersort/Qianhui_Yu_cibersort_signatures.xls")
cibersort_signatures <- cibersort_signatures[,-2]
colnames(cibersort_signatures) <- sub("fpkm_","",colnames(cibersort_signatures))
cibersort_signatures <- cibersort_signatures[order(cibersort_signatures$EnsemblID),]

dat <- cbind(EnsemblID=rownames(rpkm), as.data.frame(rpkm))
dat$EnsemblID <- as.character(dat$EnsemblID)
dat <- dat[order(dat$EnsemblID),]
dat <- dat[which(dat$EnsemblID %in% cibersort_signatures$EnsemblID),]
dat <- dat[which(dat$EnsemblID %in% rownames(cpm)),]

cibersort_signatures <- cibersort_signatures[which(cibersort_signatures$EnsemblID %in% dat$EnsemblID),]

write.table(cibersort_signatures,file="cibersort/SCZ_cibersort_signatures_fpkm.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(dat,file="cibersort/SCZ_cibersort_mixture_rpkm.txt",sep="\t",col=T,row=F,quote=F)

#- Cibersort
#- lee.marshall@vai.org
#- 100 permutations, disable quantile normalization

library(readr)
cibersort_output <- read_csv("cibersort/SCZ_cibersort_output_fpkm.csv")
covar_ciber <- cbind(covar,cibersort_output)

#- Wiskerplot
ggplot(covar_ciber, aes(DISEASE, OPC)) + 
	geom_boxplot() +
    labs(title="SCZ RNAseq STAR RSEM CIBERSORT OPC") +
	theme_bw() +
	theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
		plot.title=element_text(hjust=0.5,size=14,family="Helvetica"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		panel.background=element_blank(),
		panel.border=element_blank(),        		
		axis.title.x=element_text(size=14,family="Helvetica"),
		axis.title.y=element_text(size=14,family="Helvetica"),        
		axis.text.x=element_text(colour="black",size=10,family="Helvetica"),
		axis.text.y=element_text(colour="black",size=10,family="Helvetica"),
		axis.line=element_line(colour="black"))

astrocytes
endothelial
fetal_quiescent
microglia
neurons
oligodendrocytes
OPC


#--------------------------------- EdgeR & CIBERSORT
design1 <- model.matrix(~ GROUP + SEX + AGE + PMI + neurons, data=covar_ciber) 
#- Estimate Dispersion
y1 <- estimateDisp(y,design1)
#- quasi-likelihood F-tests
fit1 <- glmQLFit(y1,design1)
fit1 <- glmQLFTest(fit1,coef=2)
#- quasi-likelihood F-tests & FDR
qlf1 <- fit1$table
qlf1$FDR <- p.adjust(qlf1$PValue,method="BH")
#- Change Names to Entrez
#http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
qlf1$symbol <- mapIds(org.Hs.eg.db,keys=row.names(qlf1),column="SYMBOL",
	keytype="ENSEMBL",multiVals="first")
qlf1$entrez <- mapIds(org.Hs.eg.db,keys=row.names(qlf1),column="ENTREZID",
	keytype="ENSEMBL",multiVals="first")
qlf1$name <- mapIds(org.Hs.eg.db,keys=row.names(qlf1),column="GENENAME",
	keytype="ENSEMBL",multiVals="first")
#- combine gene names & cpm counts
qlf1 <- cbind(qlf1,cpm)
write.csv(qlf1, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro.csv")
#- Filter FDR <= 0.1
qlf1_FDR <- qlf1[qlf1$FDR <= 0.1,]
write.csv(qlf1_FDR, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_0.1FDR.csv")
qlf1_FDR_logFC <- qlf1_FDR[qlf1_FDR$logFC >= 1 | qlf1_FDR$logFC <= -1,]
write.csv(qlf1_FDR_logFC, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_0.1FDR_1logFC.csv")




#--------------------------------- EdgeR & CIBERSORT
design2 <- model.matrix(~ GROUP + SEX + AGE + PMI + neurons + microglia, data=covar_ciber) 
#- Estimate Dispersion
y2 <- estimateDisp(y,design2)
#- quasi-likelihood F-tests
fit2 <- glmQLFit(y2,design2)
fit2 <- glmQLFTest(fit2,coef=2)
#- quasi-likelihood F-tests & FDR
qlf2 <- fit2$table
qlf2$FDR <- p.adjust(qlf2$PValue,method="BH")
#- Change Names to Entrez
#http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
qlf2$symbol <- mapIds(org.Hs.eg.db,keys=row.names(qlf2),column="SYMBOL",
	keytype="ENSEMBL",multiVals="first")
qlf2$entrez <- mapIds(org.Hs.eg.db,keys=row.names(qlf2),column="ENTREZID",
	keytype="ENSEMBL",multiVals="first")
qlf2$name <- mapIds(org.Hs.eg.db,keys=row.names(qlf2),column="GENENAME",
	keytype="ENSEMBL",multiVals="first")
#- combine gene names & cpm counts
qlf2 <- cbind(qlf2,cpm)
write.csv(qlf2, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_micro.csv")
#- Filter FDR <= 0.1
qlf2_FDR <- qlf2[qlf2$FDR <= 0.1,]
write.csv(qlf2_FDR, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_micro_0.1FDR.csv")
qlf2_FDR_logFC <- qlf2_FDR[qlf2_FDR$logFC >= 1 | qlf2_FDR$logFC <= -1,]
write.csv(qlf2_FDR_logFC, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_micro_0.1FDR_1logFC.csv")



#--------------------------------- EdgeR & CIBERSORT
design3 <- model.matrix(~ GROUP + SEX + AGE + PMI + neurons + astrocytes + oligodendrocytes, data=covar_ciber) 
#- Estimate Dispersion
y3 <- estimateDisp(y,design3)
#- quasi-likelihood F-tests
fit3 <- glmQLFit(y3,design3)
fit3 <- glmQLFTest(fit3,coef=2)
#- quasi-likelihood F-tests & FDR
qlf3 <- fit3$table
qlf3$FDR <- p.adjust(qlf3$PValue,method="BH")
#- Change Names to Entrez
#http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
qlf3$symbol <- mapIds(org.Hs.eg.db,keys=row.names(qlf3),column="SYMBOL",
	keytype="ENSEMBL",multiVals="first")
qlf3$entrez <- mapIds(org.Hs.eg.db,keys=row.names(qlf3),column="ENTREZID",
	keytype="ENSEMBL",multiVals="first")
qlf3$name <- mapIds(org.Hs.eg.db,keys=row.names(qlf3),column="GENENAME",
	keytype="ENSEMBL",multiVals="first")
#- combine gene names & cpm counts
qlf3 <- cbind(qlf3,cpm)
write.csv(qlf3, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_astro_oligo.csv")
#- Filter FDR <= 0.1
qlf3_FDR <- qlf3[qlf3$FDR <= 0.1,]
write.csv(qlf3_FDR, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_astro_oligo_0.1FDR.csv")
qlf3_FDR_logFC <- qlf3_FDR[qlf3_FDR$logFC >= 1 | qlf3_FDR$logFC <= -1,]
write.csv(qlf3_FDR_logFC, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_PMI_neuro_astro_oligo_0.1FDR_1logFC.csv")










#-----------------------------------------------------------------------------------------
#--------------------------- SCZ RNAseq STAR EdgeR Cibersort TWO PASS NO SCHAFFOLDS ------
#-----------------------------------------------------------------------------------------
2017.12.08

#results not showing any differentially expressed genes, therefore asked shraddha for bams

mkdir Ensembl_GRCh38_p10_Softmasked_noscaffolds

wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.*

qsubs 1 24 'gunzip SAMPLENAME.fa.gz' gunzip_fasta chr_names

for i in `ls chr* | sed 's/\.fa//'` ; do sed -i '1 s/^.*$/>'${i}'/' ${i}.fa ; done

for i in `cat names` ; do cat ${i}.fa ; done > Ensembl_GRCh38_p10_Softmasked_noscaffolds.fa

#--------- STAR INDEX
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/star/2.5.3a/bin/Linux_x86_64
STAR --help
mkdir star_index_75bp

qsubs 20 24 '\
STAR --runThreadN 20 \
	--runMode genomeGenerate \
	--genomeDir star_index_75bp \
	--genomeFastaFiles Ensembl_GRCh38_p10_Softmasked_noscaffolds.fa \
	--sjdbGTFfile gencode.v26.annotation.gtf \
	--sjdbOverhang 75 \
' star_index_75bp


#--------- RSEM INDEX
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/rsem/RSEM/bin
rsem-prepare-reference --help
mkdir rsem_index

qsubs 15 24 '\
rsem-prepare-reference --num-threads 15 \
	--gtf gencode.v26.annotation.gtf  \
	Ensembl_GRCh38_p10_Softmasked_noscaffolds.fa \
	rsem_index/rsem \
' rsem_index


2017.12.11
#--------- STAR ALIGN
mkdir -p star_twopass/star
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/star/2.5.3a/bin/Linux_x86_64
STAR --help

for i in `cat names` ; do 
qsubs 4 24 'STAR --runMode alignReads \
	--runThreadN 4 \
	--genomeDir ~/labrie-primary/Bioinformatics_core/genomes/Human/Ensembl_GRCh38_p10_Softmasked_noscaffolds/star_index_75bp \
	--genomeLoad NoSharedMemory \
	--readFilesIn trimgalore/'${i}'_R1_val_1.fq.gz trimgalore/'${i}'_R2_val_2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix star_twopass/star/'${i}'_ \
	--outReadsUnmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outFilterMismatchNoverLmax 0.1 \
	--quantMode TranscriptomeSAM GeneCounts \
	--twopassMode Basic' star_aligning
done

for i in `cat names` ; do
qsubs 1 24 'samtools index star_twopass/star/'${i}'_Aligned.sortedByCoord.out.bam' samtools_index
done


#--------- RSEM COUNTING
mkdir rsem
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/rsem/RSEM/bin
rsem-calculate-expression --help

for i in `cat names` ; do 
qsubs 4 24 '\
rsem-calculate-expression \
	--paired-end \
	--num-threads 4 \
	--bam \
	--no-bam-output \
	star_twopass/star/'${i}'_Aligned.toTranscriptome.out.bam \
	~/labrie-primary/Bioinformatics_core/genomes/Human/Ensembl_GRCh38_p10_Softmasked_noscaffolds/rsem_index/rsem \
	star_twopass/rsem/'${i}'' rsem_counts
done


#-----------------------------
#----------------------------- STAR COUNTS
#-----------------------------
#--------- R EdgeR
library(readr)
covar <- read_delim("covar.tsv","\t",escape_double=FALSE,trim_ws=TRUE)

#- Merge Count data
library(edgeR)
counts <- paste(covar$ID,"_ReadsPerGene.out.tab",sep="")
counts <- readDGE(counts,path="../star",columns=c(1,4),labels=covar$ID,header=TRUE)
counts <- counts$counts[-1:-3,]
rownames(counts) <- gsub("\\..*","",rownames(counts))
write.csv(counts,file="SCZ_RNAseq_STAR_counts.csv")

#- DGEList data class
y <- DGEList(counts, group = covar$GROUP)

#- Filtering
keep <- rowSums(cpm(y)>1) >= 34
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)		#[1] 15027    34

#- Normalization
y <- calcNormFactors(y,method="TMM")
y$samples

#- CPM
cpm <- cpm(y,normalized.lib.sizes=TRUE, log=FALSE)
write.csv(cpm,file="cpm.csv")
logcpm <- cpm(y,normalized.lib.sizes=TRUE, log=TRUE)
write.csv(logcpm, file="logcpm.csv")

#- Change Names to Entrez
#http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
genes <- select(org.Hs.eg.db,keys=row.names(cpm),columns=c("SYMBOL","ENTREZID"),keytype="ENSEMBL",multiVals="first")
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes

#--------- CIBERSORT
#- Sort cibersort signature from Qianhui_Yu paper
library(readxl)
cibersort_signatures <- read_excel("cibersort/Qianhui_Yu_cibersort_signatures.xls")
cibersort_signatures <- cibersort_signatures[,-2]
colnames(cibersort_signatures) <- sub("fpkm_","",colnames(cibersort_signatures))
cibersort_signatures <- cibersort_signatures[order(cibersort_signatures$EnsemblID),]

cibersort_mixture <- cbind(EnsemblID=rownames(cpm), as.data.frame(cpm))
cibersort_mixture$EnsemblID <- as.character(cibersort_mixture$EnsemblID)
cibersort_mixture <- cibersort_mixture[order(cibersort_mixture$EnsemblID),]

#cibersort_mixture <- cibersort_mixture[which(cibersort_mixture$EnsemblID %in% cibersort_signatures$EnsemblID),]
#cibersort_signatures <- cibersort_signatures[which(cibersort_signatures$EnsemblID %in% cibersort_mixture$EnsemblID),]

write.table(cibersort_signatures,file="cibersort/SCZ_cibersort_signatures.txt",sep="\t",col=T,row=F,quote=F)
write.table(cibersort_mixture,file="cibersort/SCZ_cibersort_mixture.txt",sep="\t",col=T,row=F,quote=F)

#- Cibersort
#- lee.marshall@vai.org
#- 100 permutations, disable quantile normalization

library(readr)
cibersort_output <- read_csv("cibersort/SCZ_cibersort_output.csv")
covar_ciber <- cbind(covar,cibersort_output)

#- Wiskerplot
library(ggplot2)

ggplot(covar_ciber, aes(DISEASE,astrocytes)) + 
	geom_boxplot() +
    labs(title="SCZ RNAseq STAR RSEM CIBERSORT") +
	theme_bw() +
	theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
		plot.title=element_text(hjust=0.5,size=14,family="Helvetica"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		panel.background=element_blank(),
		panel.border=element_blank(),        		
		axis.title.x=element_text(size=14,family="Helvetica"),
		axis.title.y=element_text(size=14,family="Helvetica"),        
		axis.text.x=element_text(colour="black",size=10,family="Helvetica"),
		axis.text.y=element_text(colour="black",size=10,family="Helvetica"),
		axis.line=element_line(colour="black"))


astrocytes
endothelial
fetal_quiescent
microglia
neurons
oligodendrocytes
OPC

#----------------------------- EdgeR & CIBERSORT
design1 <- model.matrix(~ GROUP + SEX + AGE + neurons, data=covar_ciber) 
#- Estimate Dispersion
y1 <- estimateDisp(y,design1)
#- generalized linear model
fit1 <- glmFit(y1,design1)
#- liklihood ratio test
lrt1 <- glmLRT(fit1,coef=2)
lrt1 <- lrt1$table
lrt1$FDR <- p.adjust(lrt1$PValue,method="BH")
#- combine gene names & cpm counts
lrt1 <- cbind(lrt1,fit1$genes,cpm)
write.csv(lrt1, file="SCZ_RNAseq_STAR_EDGER_glmlrt_Cond_Sex_Age_neuro.csv")
#- Filter FDR <= 0.1
lrt1_FDR <- lrt1[lrt1$FDR <= 0.1,]
write.csv(lrt1_FDR, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_neuro_0.1FDR.csv")

design2 <- model.matrix(~ GROUP + SEX + AGE + neurons + astrocytes + oligodendrocytes + OPC, data=covar_ciber) 
#- Estimate Dispersion
y2 <- estimateDisp(y,design2)
#- generalized linear model 
fit2 <- glmFit(y2,design2)
#- liklihood ratio test
lrt2 <- glmLRT(fit2,coef=2)
lrt2 <- lrt2$table
lrt2$FDR <- p.adjust(lrt2$PValue,method="BH")
#- combine gene names & cpm counts
lrt2 <- cbind(lrt2,fit2$genes,cpm)
write.csv(lrt2, file="SCZ_RNAseq_STAR_EDGER_glmlrt_Cond_Sex_Age_neuro_astro_oligo_opc.csv")
#- Filter FDR <= 0.1
lrt2_FDR <- lrt2[lrt1$FDR <= 0.1,]
write.csv(lrt2_FDR, file="SCZ_RNAseq_STAR_RSEM_EDGER_Cond_Sex_Age_neuro_astro_oligo_opc_0.1FDR.csv")




#----------------------------- LIMMA VOOM on STAR COUNTS
#--------- Design
design1 <- model.matrix(~GROUP + SEX + AGE + neurons, data=covar_ciber) 
colnames(design1) <- gsub("GROUP","",colnames(design1))
#---------  Voom
v1 <- voom(y,design1,plot=TRUE)
#---------  Limma
vfit <- lmFit(v1,design1)
efit <- eBayes(vfit)
limmavoom <- toptable(efit,coef=2,adjust.method="BH",sort="none",n=Inf)
limmavoom <- cbind(limmavoom,efit$genes,cpm)
write.csv(limmavoom, file="SCZ_RNAseq_STAR_limma_voom_Cond_Sex_Age_neuro.csv")
limmavoom_FDR <- limmavoom[limmavoom$adj.P.Val <= 0.1,]
write.csv(limmavoom_FDR, file="SCZ_RNAseq_STAR_limma_voom_Cond_Sex_Age_neuro_0.1FDR.csv")

#--------- Design
design2 <- model.matrix(~GROUP + SEX + AGE + neurons + astrocytes + oligodendrocytes + OPC, data=covar_ciber) 
colnames(design2) <- gsub("GROUP","",colnames(design2))
#---------  Voom
v2 <- voom(y,design2,plot=TRUE)
#---------  Limma
vfit2 <- lmFit(v2,design2)
efit2 <- eBayes(vfit2)
limmavoom2 <- toptable(efit2,coef=2,adjust.method="BH",sort="none",n=Inf)
limmavoom2 <- cbind(limmavoom2,efit2$genes,cpm)
write.csv(limmavoom2, file="SCZ_RNAseq_STAR_limma_voom_Cond_Sex_Age_neuro_astro_oligo_opc.csv")
limmavoom2_FDR <- limmavoom2[limmavoom2$adj.P.Val <= 0.1,]
write.csv(limmavoom2_FDR, file="SCZ_RNAseq_STAR_limma_voom_Cond_Sex_Age_neuro_astro_oligo_opc_0.1FDR.csv")




#-----------------------------
#----------------------------- RSEM COUNTS
#-----------------------------
#--------- R EdgeR
library(readr)
covar <- read_delim("covar.tsv","\t",escape_double=FALSE,trim_ws=TRUE)

#- Merge Count data
library(edgeR)
counts <- paste(covar$ID,".genes.results",sep="")
counts <- readDGE(counts,path="../rsem",columns=c(1,5),labels=covar$ID,header=TRUE)
counts <- counts$counts
rownames(counts) <- gsub("\\..*","",rownames(counts))
write.csv(counts,file="SCZ_RNAseq_RSEM_counts.csv")

fpkm_counts <- paste(covar$ID,".genes.results",sep="")
fpkm_counts <- readDGE(fpkm_counts,path="../rsem",columns=c(1,7),labels=covar$ID,header=TRUE)
fpkm_counts <- fpkm_counts$counts
rownames(fpkm_counts) <- gsub("\\..*","",rownames(fpkm_counts))
write.csv(fpkm_counts,file="SCZ_RNAseq_RSEM_fpkm_counts.csv")

#- DGEList data class
y <- DGEList(counts, group = covar$GROUP)

#- Filtering
keep <- rowSums(cpm(y)>1) >= 34
y <- y[keep, ,keep.lib.sizes=FALSE]
dim(y)		#[1] 15308    34

#- Normalization
y <- calcNormFactors(y,method="TMM")
y$samples

#- CPM
cpm <- cpm(y,normalized.lib.sizes=TRUE, log=FALSE)
write.csv(cpm,file="SCZ_RNAseq_RSEM_cpm.csv")
logcpm <- cpm(y,normalized.lib.sizes=TRUE, log=TRUE)
write.csv(logcpm, file="SCZ_RNAseq_RSEM_logcpm.csv")

#- Change Names to Entrez
#http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
genes <- select(org.Hs.eg.db,keys=row.names(cpm),columns=c("SYMBOL","ENTREZID"),keytype="ENSEMBL",multiVals="first")
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes

#--------- CIBERSORT
#- Sort cibersort signature from Qianhui_Yu paper
library(readxl)
cibersort_signatures <- read_excel("cibersort/Qianhui_Yu_cibersort_signatures.xls")
cibersort_signatures <- cibersort_signatures[,-2]
colnames(cibersort_signatures) <- sub("fpkm_","",colnames(cibersort_signatures))
cibersort_signatures <- cibersort_signatures[order(cibersort_signatures$EnsemblID),]

cibersort_mixture <- cbind(EnsemblID=rownames(fpkm_counts), as.data.frame(fpkm_counts))
cibersort_mixture <- cibersort_mixture[keep,]
cibersort_mixture$EnsemblID <- as.character(cibersort_mixture$EnsemblID)
cibersort_mixture <- cibersort_mixture[order(cibersort_mixture$EnsemblID),]

write.table(cibersort_signatures,file="cibersort/SCZ_cibersort_signatures_fpkm.txt",sep="\t",col=T,row=F,quote=F)
write.table(cibersort_mixture,file="cibersort/SCZ_cibersort_mixture_fpkm.txt",sep="\t",col=T,row=F,quote=F)

#- Cibersort
#- lee.marshall@vai.org
#- 100 permutations, disable quantile normalization

library(readr)
cibersort_output <- read_csv("cibersort/SCZ_cibersort_output_fpkm.csv")
covar_ciber <- cbind(covar,cibersort_output)


#----------------------------- EdgeR & CIBERSORT on RSEM COUNTS
design1 <- model.matrix(~ GROUP + SEX + AGE + neurons, data=covar_ciber) 
#- Estimate Dispersion
y1 <- estimateDisp(y,design1)
#- generalized linear model
fit1 <- glmFit(y1,design1)
#- liklihood ratio test
lrt1 <- glmLRT(fit1,coef=2)
lrt1 <- lrt1$table
lrt1$FDR <- p.adjust(lrt1$PValue,method="BH")
#- combine gene names & cpm counts
lrt1 <- cbind(lrt1,fit1$genes,cpm)
write.csv(lrt1, file="SCZ_RNAseq_RSEM_EDGER_glmlrt_Cond_Sex_Age_neuro.csv")
#- Filter FDR <= 0.1
lrt1_FDR <- lrt1[lrt1$FDR <= 0.1,]
write.csv(lrt1_FDR, file="SCZ_RNAseq_RSEM_EDGER_glmlrt_Cond_Sex_Age_neuro_0.1FDR.csv")

design2 <- model.matrix(~ GROUP + SEX + AGE + neurons + astrocytes + oligodendrocytes, data=covar_ciber) 
#- Estimate Dispersion
y2 <- estimateDisp(y,design2)
#- generalized linear model 
fit2 <- glmFit(y2,design2)
#- liklihood ratio test
lrt2 <- glmLRT(fit2,coef=2)
lrt2 <- lrt2$table
lrt2$FDR <- p.adjust(lrt2$PValue,method="BH")
#- combine gene names & cpm counts
lrt2 <- cbind(lrt2,fit2$genes,cpm)
write.csv(lrt2, file="SCZ_RNAseq_RSEM_EDGER_glmlrt_Cond_Sex_Age_neuro_astro_oligo.csv")
#- Filter FDR <= 0.1
lrt2_FDR <- lrt2[lrt1$FDR <= 0.1,]
write.csv(lrt2_FDR, file="SCZ_RNAseq_RSEM_EDGER_glmlrt_Cond_Sex_Age_neuro_astro_oligo_0.1FDR.csv")


#----------------------------- LIMMA VOOM on RSEM COUNTS
#--------- Design
design1 <- model.matrix(~GROUP + SEX + AGE + neurons, data=covar_ciber) 
colnames(design1) <- gsub("GROUP","",colnames(design1))
#---------  Voom
v1 <- voom(y,design1,plot=TRUE)
#---------  Limma
vfit <- lmFit(v1,design1)
efit <- eBayes(vfit)
limmavoom <- toptable(efit,coef=2,adjust.method="BH",sort="none",n=Inf)
limmavoom <- cbind(limmavoom,efit$genes,cpm)
write.csv(limmavoom, file="SCZ_RNAseq_RSEM_limma_voom_Cond_Sex_Age_neuro.csv")
limmavoom_FDR <- limmavoom[limmavoom$adj.P.Val <= 0.1,]
write.csv(limmavoom_FDR, file="SCZ_RNAseq_RSEM_limma_voom_Cond_Sex_Age_neuro_0.1FDR.csv")

#--------- Design
design2 <- model.matrix(~GROUP + SEX + AGE + neurons + astrocytes + oligodendrocytes, data=covar_ciber) 
colnames(design2) <- gsub("GROUP","",colnames(design2))
#---------  Voom
v2 <- voom(y,design2,plot=TRUE)
#---------  Limma
vfit2 <- lmFit(v2,design2)
efit2 <- eBayes(vfit2)
limmavoom2 <- toptable(efit2,coef=2,adjust.method="BH",sort="none",n=Inf)
limmavoom2 <- cbind(limmavoom2,efit2$genes,cpm)
write.csv(limmavoom2, file="SCZ_RNAseq_RSEM_limma_voom_Cond_Sex_Age_neuro_astro_oligo.csv")
limmavoom2_FDR <- limmavoom2[limmavoom2$adj.P.Val <= 0.1,]
write.csv(limmavoom2_FDR, file="SCZ_RNAseq_RSEM_limma_voom_Cond_Sex_Age_neuro_astro_oligo_0.1FDR.csv")





#-----------------------------
#----------------------------- SALMON COUNTS
#-----------------------------
#--------- Salmon Quant
export PATH=$PATH:~/labrie-primary/Bioinformatics_core/software/salmon/0.8.2/bin/
salmon quant --help-reads
bootstrapping for genes, gibbs sampling for isoforms 

for i in `cat names` ; do
qsubs 4 24 "salmon quant \
	--index ~/labrie-primary/Bioinformatics_core/transcriptomes/GRCh38_p10_gencode_v26/Salmon_index \
	--libType A \
	--mates1 trimgalore/${i}_R1_val_1.fq.gz \
	--mates2 trimgalore/${i}_R2_val_2.fq.gz \
	--output salmon_quasi/${i} \
	--threads 4 \
	--geneMap ~/labrie-primary/Bioinformatics_core/transcriptomes/GRCh38_p10_gencode_v26/GRCh38_p10_gencode_v26.gtf \
	--useVBOpt \
	--numBootstraps 100" "salmon_quasi"
done





