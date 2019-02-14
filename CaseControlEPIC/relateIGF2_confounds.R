# relate IGF2 methylation to lifestyle variables (smoking, medication, brain wt)
rm(list=ls())

require(ggplot2)

inFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/preprocessing/CaseControlEPIC_CLEAN_171127.Rdata"
outDir <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/output_files/CaseControlEPIC/SUS19398/dmp_DX_QTL/agesexPMI_PC12_171129"
phenoFile <- "/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/input_files/NARSAD_sampleKey_allMeds_180103.txt"

p2get <- c("cg07096953","cg02613624","cg22956483")

load(inFile)

require(minfi)

# integrate with confounds
metadat <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
colnames(metadat)[2] <- "Sample_ID"
colnames(metadat)[which(colnames(metadat)=="TissueBank.ID")] <- "ID"
colnames(metadat)[which(colnames(metadat)=="Smoker.status")] <- "smoker.status"
colnames(metadat)[which(colnames(metadat)=="Parsed.antipsychotics")] <- "antipsych.parsed"
colnames(metadat)[which(colnames(metadat)=="mood.stabilizer..valproate")] <- "valproate"
colnames(metadat)[which(colnames(metadat)=="mood.stabilizer..lithium")] <- "lithium"

idx <- which(metadat$valproate == "yes")
metadat$antipsych.parsed[idx] <- sprintf("%s;valproate", 
	metadat$antipsych.parsed[idx])
idx <- which(metadat$lithium == "yes")
metadat$antipsych.parsed[idx] <- sprintf("%s;lithium", 
	metadat$antipsych.parsed[idx])

metadat <- metadat[,c("ID","Sample_ID","BrainWt","smoker.status","antipsych.parsed","valproate","lithium")]

# combine methylation with patient ID
pd <- pData(MSet.genome)
idx <- !duplicated(pd$ID)

MSet.genome <- MSet.genome[,idx]
pd <- pData(MSet.genome)
pd$Sample_Rowname <- rownames(pd)
merged <- merge(x=pd,y=metadat,by="Sample_ID",all.x=TRUE)
midx <- match(rownames(pd),merged$Sample_Rowname)
if (all.equal(merged$Sample_Rowname[midx],rownames(pd))!=TRUE) {
	cat("ids don't match\n"); browser()
}

if (nrow(pd)!=nrow(merged)) {
	cat("Some methylation samples not found in pheno table\n");
	browser()
}

pd <- merged[midx,]

betas <- as.data.frame(t(getBeta(MSet.genome[p2get,])))

if (all.equal(rownames(betas),pd$Sample_Rowname)!=TRUE) {
	cat("Betas and pd don't match"); browser()
}

betas$ID <- pd$ID.y
betas$DX <- pd$DX
betas$AGE <- pd$AGE
betas$SEX <- pd$SEX
betas$PMI <- pd$PMI
# add metadata columns
betas$BrainWt <- pd$BrainWt

pd$valproate[which(pd$valproate!="yes")] <- "unlisted"
betas$valproate <- pd$valproate
pd$lithium[which(pd$lithium!="yes")] <- "unlisted"
betas$lithium <- pd$lithium

# binarize antipsychotic use
apsych <- pd$antipsych.parsed
apsych[which(apsych!="none")] <- "some"
betas$antipsych <- apsych
# binarize smoking status
smk <- rep(NA,nrow(betas))
pd$smoker.status <- trimws(pd$smoker.status)
smk[which(pd$smoker.status %in% c("non smoker","No","nonsmoker"))] <- "Nonsmoker"
smk[which(pd$smoker.status %in% c("Yes","previous smoker, none at tod","smoker"))] <- "Smoker"
betas$smoker.status <- smk

betas$C1 <- pd$C1
betas$C2 <- pd$C2

# remove duplicates
cat(sprintf("Removing duplicate sample IDs\n"))
betas <- betas[!duplicated(betas$ID),]

betas$antipsych <- factor(betas$antipsych,levels=c("none","some"))
betas$smoker.status <- factor(betas$smoker.status,
	levels=c("Nonsmoker","Smoker"))
betas$DX <- factor(betas$DX,levels=c("control","case"))
betas$PMI <- as.numeric(betas$PMI)
# plot confounds for selected probe
plotMvsConfound <- function(probeName) {
	p <- list()
	betas2 <- betas[,c(probeName,"DX","SEX","AGE","BrainWt",
		"smoker.status","antipsych","PMI","valproate","lithium","C1","C2")]
	colnames(betas2)[1] <- "probe"

	cat("---------------------------\n")
	cat(sprintf("%s: effect of DX, control smoking;antipsych;PMI\n",
		probeName))
	m1 <- lm(probe~antipsych+smoker.status+AGE+SEX+PMI+C1+C2,data=betas2)
	m2 <- lm(probe~DX+antipsych+smoker.status+AGE+SEX+PMI+C1+C2,data=betas2)
	res <- anova(m1,m2)
	pval <- res[["Pr(>F)"]]
	cat(sprintf("%s: pvalue from nested AOV = %1.2e\n", probeName,pval[2]))
	cat("---------------------------\n")

	# brain weight
	p[[1]] <- ggplot(betas2,aes(y=probe,x=BrainWt)) 
	p[[1]] <- p[[1]] + geom_point(aes(colour=DX,pch=SEX)) + geom_smooth(method=lm)
	
	# smoker status
	cat("* Smoker status (ever/never)\n")
	p[[2]] <- ggplot(betas2,aes(smoker.status,probe))
	p[[2]] <- p[[2]]+geom_boxplot(aes(colour=DX))
	p[[2]] <- p[[2]] + ggtitle("probe:Smoking status")
	
	#  antipsychotics taken
	cat("* Antipsych status (ever/never)\n")
	p[[3]] <- ggplot(betas2,aes(antipsych,probe))+geom_boxplot(aes(colour=DX))
	p[[3]] <- p[[3]] + ggtitle("probe: Antipsychotics")
	cat(sprintf("Sample breakdown for NONE\n"))
	
	# pmi
	cat("* PMI\n")
	p[[4]] <- ggplot(betas2,aes(PMI,probe))+geom_point(aes(colour=DX))+geom_smooth(aes(colour=DX),method="lm")
	p[[4]] <- p[[4]] + ggtitle("probe: PMI")

	# valproate
	p[[5]] <- ggplot(betas2,aes(valproate,probe))+geom_boxplot(aes(colour=DX))
	p[[5]] <- p[[5]] + ggtitle("Valproate")

	# lithium
	p[[6]] <- ggplot(betas2,aes(lithium,probe))+geom_boxplot(aes(colour=DX))
	p[[6]] <- p[[6]] + ggtitle("Lithium")

	return(list(plotList=p,res=pval))
}


source("multiplot.R")
dt <- format(Sys.Date(),"%y%m%d")
logFile <- sprintf("%s/igf2_confound_%s.log",outDir,dt)
sink(logFile,split=TRUE)
tryCatch({
for (probeName in p2get) {
	print(probeName)
	pdfFile <- sprintf("%s/igf2_confound_%s_%s.pdf",outDir,probeName,dt)
	res <- plotMvsConfound(probeName)
	pdf(pdfFile,width=11,height=6)
	multiplot(plotlist=res[["plotList"]],layout=matrix(1:6,ncol=3))
	dev.off()
}
cat("Case/control breakdown for no antipsychotic\n")
	print(table(betas$DX[which(betas$antipsych=="none")],useNA="always"))
cat("Breakdown for smoking status\n")
	print(table(betas$smoker.status),useNA="always")
cat("Breakdown for valproate\n")
	print(table(betas$valproate),useNA="always")
cat("Breakdown for lithium\n")
	print(table(betas$lithium),useNA="always")

cat("Test against average methylation\n")
betas$M <- rowMeans(betas[,1:3])
	m1 <- lm(M~antipsych+smoker.status+AGE+SEX+PMI+C1+C2,data=betas)
	m2 <- lm(M~DX+antipsych+smoker.status+AGE+SEX+PMI+C1+C2,data=betas)
	res <- anova(m1,m2)
	pval <- res[["Pr(>F)"]]
	cat(sprintf("avg: pvalue from nested AOV = %1.2e\n", pval[2]))
},error=function(ex){print(ex)},finally={sink(NULL)})

