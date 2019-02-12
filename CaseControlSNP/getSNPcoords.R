#' get coords for pgc snps

pgcFile<-"/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/pgc/scz2_plt1e-9.txt"

require(SNPlocs.Hsapiens.dbSNP144.GRCh37)
dat <- read.delim(pgcFile,sep="\t",h=F,as.is=T)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
rsid <- dat[grep("^rs",dat[,1]),1]
cat(sprintf("%i snps\n",length(rsid)))

cat("* Fetching coords\n")
 df <- cbind(rsid=rsid[1:10],
			 chrom=names(locs), 
			 cstart=locs)


# get coordinates in parallel chunks of 100s.
# SP ran timing and running snpid2loc() call one at a time
# is faster
require(foreach)
require(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

t0 <- Sys.time();
out <- foreach(spos=seq(200,300,100),
			   .packages="SNPlocs.Hsapiens.dbSNP144.GRCh37") %dopar% {
	chrom <- character(100)
	coord <- integer(100)
	epos <- spos+99
	if (epos > length(rsid)) epos <- length(rsid)
	ctr<-1
	for (k in spos:epos) { 
		x <- snpid2loc(snps,rsid[k]); 
		coord[ctr] <- x; 
		chrom[ctr] <- names(x)
		ctr <- ctr+1
	} 
	x <- cbind(rsid[spos:epos],names(coord),coord)
	x
}
print(Sys.time()-t0)
out2 <- do.call("rbind", out)

