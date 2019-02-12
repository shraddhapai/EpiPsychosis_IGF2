# run plink on genotype data to compute lambda
# use age and pmi as covariates

plink <- "~/software/plink/plink"

inFile <- "/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/SUS19399_New/hg19/SUS19399.hg19.sorted_-clean"
phenoFile <- "/home/shraddhapai/Epigenetics/NARSAD/input_files/NARSAD_sampleKey_171129.txt"
pcFile <-sprintf("%s_PCA.eigenvec",inFile) 
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/computeLambda_%s",dirname(inFile),dt)

if (!file.exists(outDir)) dir.create(outDir,recursive=FALSE)

# remove tech reps
fam <- read.table(sprintf("%s.fam",inFile),sep=" ",h=F,as.is=T)
idx <- grep("-R[23]",fam[,2])
cat(sprintf("removing %i tech reps\n",length(idx)))
fam <- fam[-idx,]

# now match with pheno
pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
pheno <- pheno[,c("Internal.ID","AGE","SEX","PMI","DIST.DX")]
dx <- rep("Control",nrow(pheno))
dx[which(pheno$DIST.DX %in% c("Bipolar","Schizophrenia"))] <- "Case"
pheno$DX <- dx

id <- sub("-R1","",fam[,2])
tokeep <- intersect(id,pheno$Internal.ID)
fam <- fam[which(id %in% tokeep),1:2]
cat(sprintf("Keeping %i samples intersecting with pheno\n",nrow(fam)))

keepFile <- sprintf("%s/tokeep.txt",outDir)
write.table(fam,file=keepFile,sep=" ",col=F,row=F,quote=F)
genoFile <- sprintf("%s/inputSamps",outDir)
cmd <- sprintf("%s --bfile %s --keep %s --make-bed --out %s",
	plink, inFile, keepFile,genoFile)
cat("removing tech reps from plink\n")
system(cmd)

# set up covariate file and diagnosis
fam <- read.table(sprintf("%s/inputSamps.fam",outDir),sep=" ",h=F,as.is=T)
pcdat <- read.delim(pcFile,sep=" ",h=F,as.is=T)
midx <- match(fam[,2],pcdat[,2])
if (all.equal(pcdat[midx,2],fam[,2])!=TRUE) {
	cat("pca and fam don't match"); browser()
}
pcdat <- pcdat[midx,]

id <- sub("-R1","",fam[,2])
pheno[,1] <- as.character(pheno[,1])
midx <- match(id,pheno$Internal.ID)
if (all.equal(pheno$Internal.ID[midx],id)!=TRUE) {
	cat("fam/pheno don't match\n")
	browser()
}
pheno <- pheno[midx,]
pheno$FID <- fam[,1]
pheno$IID <- fam[,2]
sex <- rep(2,nrow(pheno))
sex[which(pheno$SEX=="Male")] <- 1
dx <- rep(2,nrow(pheno))
dx[which(pheno$DX=="Control")] <- 1

fam[,5] <- sex
fam[,6] <- dx
pheno$FID <- fam[,1]
pheno$IID <- fam[,2]
pheno$C1 <- pcdat[,3]
pheno$C2 <- pcdat[,4]

covDat <- pheno[,c("FID","IID","AGE","PMI","C1","C2")]
write.table(fam,file=sprintf("%s.fam",genoFile),sep=" ",col=F,row=F,quote=F)
write.table(covDat,file=sprintf("%s.cov",genoFile),sep="\t",col=F,row=F,quote=F)

# call gwas
cat("running assoc")
cmd <- sprintf("%s --bfile %s --fam %s.fam --logistic --covar %s.cov --adjust --out %s_binAssoc",
	plink,genoFile,genoFile,genoFile,genoFile)
system(cmd)

