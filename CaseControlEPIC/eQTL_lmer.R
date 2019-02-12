#' performs eQTL analysis using lm
#'
#' @param snpFile (char) path to .traw file containing snps for which
#' eQTL analysis is to be performed
#' @param M (matrix) beta values. Rows are cpgs, columns are samples
#' @param mLocs (GRanges) locations of CpG probes
#' @param cvrt (data.frame) rows are samples columns are covariates
#' must have AGE, SEX, DX, C1 and C2 (pc1,2 for genetic ancestry)
#' Sample ID should be in column "Sample_Name"
#' @param addGenoDiseaseTerm (logical) if TRUE adds an interaction
#' term to con
#' @param lim2snps (char) if "*" starts with all snps in snpFile. Otherwise
#' expects a character vector with ids of snps to include
#' @param MAX_DIST (integer) max distance of SNP-DMP to run meQTL analysis.
#' set to Inf for all-to-all comparison (transeQTL)
eQTL_lmer <- function(snpFile, M,mLocs,cvrt,addGenoDiseaseTerm=FALSE,
	lim2snps="*",MAX_DIST=1e6) {
	require(lme4)

	source("meQTL_prepareGeno.R")
	snpDat <- meQTL_prepareGeno(snpFile=snpFile)
	snpGR <- snpDat$snpGR
	snps <- snpDat$geno

	both <- intersect(colnames(snps), colnames(M))
	
	snps	<- snps[,which(colnames(snps) %in% both)]
	M		<- M[,which(colnames(M) %in% both),drop=FALSE]
	cvrt	<- cvrt[which(cvrt$Sample_Name %in% both),]

	cat(sprintf("\tAfter merging, %i samples left\n",length(both)))
	midx <- match(colnames(M), colnames(snps))
	if (all.equal(colnames(snps)[midx],colnames(M))!=TRUE) {
		cat("snp/xpr don't match\n"); browser()
	}
	snps <- snps[,midx]

	midx <- match(colnames(M), cvrt$Sample_Name)
	if (all.equal(colnames(M),cvrt$Sample_Name[midx])!=TRUE) {
		cat("genes and cvrt sample order doesn't match\n")
		browser()
	}
	cvrt <- cvrt[midx,]

	if (all.equal(colnames(snps),cvrt$Sample_Name)!=TRUE) {
		cat("snps and cvrt sample order doesn't match\n")
		browser()
	}

# remove technical replicates
msamp <- sub("-R[123456]","",colnames(M))
newM <- c()
for (k in 1:nrow(M)) {
	tmp <- aggregate(M[k,],by=list(Sample=msamp),FUN=mean)
	newM <- rbind(newM, tmp[,2])
	if (k>1) {
		if (all.equal(colnames(newM),tmp[,1])!=TRUE) { 
		cat("samp order changed.\n"); browser()
		}
	} else {
		colnames(newM) <- tmp[,1]
	}
}

rownames(newM) <- rownames(M); 
oldM <- M
M <- newM; rm(newM)
# for genotypes, take first one
idx <- grep("-R[23456]",colnames(snps))
snps <- snps[,-idx]
colnames(snps) <- sub("-R1","",colnames(snps))

# now align
M <- M[,order(colnames(M)),drop=FALSE]
snps <- snps[,order(colnames(snps))]
cvrt <- as.data.frame(cvrt)
cvrt <- cvrt[-grep("-R[23456]",cvrt$Sample_Name),]
cvrt$Sample_Name <- sub("-R1","",cvrt$Sample_Name)
cvrt <- cvrt[order(cvrt$Sample_Name),]

if (all.equal(cvrt$Sample_Name, colnames(snps))!=TRUE) {
	cat("samp names for cvrt and snps don't match\n");browser()
}
if (all.equal(cvrt$Sample_Name, colnames(M))!=TRUE) {
	cat("samp names for cvrt and M don't match\n"); browser();
}

# old -- ignoring tech reps now
###	idx <- grep("-R[23456]",cvrt$Sample_Name)
###	cat(sprintf("Got %i tech reps\n",length(idx)))
###	sampName <- sub("-R[123456]","",cvrt$Sample_Name)
###	cvrt$UqSamp <- factor(sampName)
###

	cvrt$SEX <- factor(cvrt$SEX)
	cvrt$DX <- factor(cvrt$DX,levels=c("control","case"))
	
	cvrt <- cvrt[,c("Sample_Name","DX","SEX","AGE","C1","C2","PMI")]

	cat(sprintf("Processing %i SNPs\n", nrow(snps)))
	cat(sprintf("Have %i samples\n", ncol(snps)))

	if (lim2snps[1]!="*") {
		cat("\tSNP filter provided!") 
		snps <- snps[which(rownames(snps)%in% lim2snps),]
		snpGR <- snpGR[which(snpGR$name %in% lim2snps)]
		cat(sprintf(" Limiting to %i SNPs\n", length(snpGR)))
	}
	ol <- NULL
	if (!is.infinite(MAX_DIST)) {
		zoneGR <- resize(mLocs,fix="center",width=MAX_DIST)
		ol <- findOverlaps(snpGR,zoneGR)
		idx <-unique(queryHits(ol))
		snps <- snps[idx,]
		snpGR <- snpGR[idx]
		ol <- findOverlaps(snpGR,zoneGR)
		cat(sprintf("dist filter: %i left\n", nrow(snps)))
	}

	# remove SNPs that don't have min num samp per genotype
	s1 <- rowSums(snps==0)
	s2 <- rowSums(snps==1)
	s3 <- rowSums(snps==2)
	idx <- which(s1>=10 & s2>=10 & s3>=10)
	snpGR <- snpGR[idx]
	snps <- snps[idx,]
	cat(sprintf("geno count filter: %i left\n", length(snpGR)))

cvrt_main <- cvrt; 
rm(cvrt)
res <- list()
res_ctr <- 1
require(parallel)
require(foreach)
require(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)

t0 <-Sys.time()
res <- foreach (snp_idx=1:nrow(snps),.packages=c("lme4","GenomicRanges")) %dopar% {
	cvrt <- cvrt_main;
	cvrt$GENO <- as.integer(snps[snp_idx,])
	snp_name <- rownames(snps)[snp_idx]

	if (!is.infinite(MAX_DIST)) {
		myM <- subsetByOverlaps(zoneGR,snpGR[snp_idx])
	} else { # don't filter by distance
		myM <- mLocs
	}
	curM <- names(myM);rm(myM)

tmp <- list(); ctr <-1
for (nm in curM) {
	print(nm)
	if (!is.infinite(MAX_DIST)) {
	d <- distanceToNearest(mLocs[nm],snpGR[snp_idx])
	d <- elementMetadata(d)@listData$distance
		if (abs(d) >= MAX_DIST/2) {
		cat("TOO FAR!\n"); browser()
		}
	}
	cvrt_tmp <- cvrt; cvrt_tmp$M <- M[nm,]
	tmp_res <- c()
	if (addGenoDiseaseTerm) {
		fit <- summary(lm(M~GENO+DX+(GENO:DX)+AGE+SEX,data=cvrt_tmp))
		tmp_res <- coef(fit)["GENO:DXcase",c(1,4)]
	} else {
		fit <- summary(lm(M~GENO+DX+AGE+SEX+C1+C2,data=cvrt_tmp))
		tmp_res <- coef(fit)["GENO",c(1,4)]
	}
	res_cur <- cbind(nm,snp_name,tmp_res[1],tmp_res[2])
	tmp[[ctr]] <-res_cur
	ctr <- ctr+1
}
tmp <- do.call("rbind",tmp)
tmp
} 
stopCluster(cl)
t1 <- Sys.time()
cat(sprintf("Loop took %1.2f s\n", t1-t0))
tmp <- res
res <- do.call("rbind",res)
colnames(res) <- c("cpg","snps","slope","p")
res <- data.frame(res)
res[,1] <- as.character(res[,1])
res[,2] <- as.character(res[,2])
res[,3] <- as.numeric(as.character(res[,3]))
res[,4] <- as.numeric(as.character(res[,4]))

comb <- res
cat("* Combining results\n")
cat(sprintf("Num tests = %i  (%i unique SNPs; %i CpGs)\n",
			nrow(comb), sum(!duplicated(comb$snps)),
			sum(!duplicated(comb$cpg))))
comb$FDR <- p.adjust(comb$p,method="BH")
bonf <- 0.05/nrow(comb)
comb$PASS_BONF <- comb$p < bonf
comb <- comb[order(comb$FDR),]

return(list(res=comb,snps=snps,M=M,cvrt=cvrt_main))
}
