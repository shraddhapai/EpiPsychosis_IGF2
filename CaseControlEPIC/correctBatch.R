# run sva to estimate batch effect

require(sva)
#' @param mset
#' @param runCombat (logical) if TRUE runs combat using array number ("Slide")
#' @return if runCombat=TRUE, returns a list with:
#' 1) combat_Mset: (GenomicRatioSet) GenomicRatioSet with combat-corrected beta values.
#' 2) f_p: (numeric) vector of pvalues for F-test (result of f.pvalue)
#' 3) f_q: (numeric) vector of Qvalues ; f_p after BH correction
correctBatch <- function(mset,runCombat=FALSE){
betas <- getBeta(mset)
pd		<- pData(mset)

pd$AGE	<- as.integer(pd$AGE)
pd$DX	<- factor(pd$DX, levels=c("case","control")) 

cat("Examining variables\n")
cat("\nAge\n")
print(summary(pd$AGE))

cat("\nDiagnosis\n")
print(table(pd$DX))

cat("\n* Computing surrogate variables\n")
mod <- model.matrix(~1+AGE+factor(SEX)+DX+PMI,data=pd)
mod0 <- model.matrix(~1+AGE+factor(SEX)+PMI,data=pd)
n.sv <- sva::num.sv(betas,mod,method="leek",vfilter=100*1000)
cat(sprintf("\tNum surrogate variables= %i\n",n.sv))

if (n.sv >=1 )  {
	cat("perform sva - not implemented\n")
	browser()
}

out <- TRUE
# run combat
if (runCombat) {
	cat("* Using ComBaT to correct for array number\n")
	modcombat <- mod0 
	# combat_edata contains batch-adjusted beta values
	combat_edata <- sva::ComBat(dat=betas,batch=pd$Slide, 
						   mod=modcombat, 
						   par.prior=TRUE, prior.plots=TRUE)

	combat_mset <- GenomicRatioSet(gr=getLocations(mset), 
								   Beta=combat_edata, 
								   pData=pData(mset), 
								   annotation=annotation(mset), 
								   preprocessMethod=preprocessMethod(mset))

	# F-test to see affect of age (mod) after accounting for null model (mod0).
	# use batch-corrected input (combat_edata).
	pValuesComBat <- sva::f.pvalue(combat_edata,mod,mod0)
	qValuesComBat <- p.adjust(pValuesComBat,method="BH")
	cat(sprintf("%i loci have Q < 0.05 (effect of diagnosis)\n",
				sum(qValuesComBat < 0.05))) 
	rm(mod)
	par(mfrow=c(1,2))
	hist(pValuesComBat,n=100,main="combat pvalues")
	hist(qValuesComBat,n=100,main="combat qvalues")

	out <- list(combat_mset=combat_mset,f_p=pValuesComBat,
				f_q=qValuesComBat)
}
	out
}
