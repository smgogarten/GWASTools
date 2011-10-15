
qualityScoreBySnp<- 
function(intenData, 
         genoData, 
         scan.exclude = NULL,
         block.size = 5000, 
         verbose = TRUE) 
{

	# block.size is the number of rows to read and process at one time (i.e. one pass of the loop below)
	# scan.exclude is a vector of integer ids of samples to exclude from the missing fraction calculation (use NULL if none to exclude)

        # check that intenData has quality
        if (!hasQuality(intenData)) stop("quality not found in intenData")
  
        # check that dimensions of intenData and genoData are equal
        intenSnpID <- getSnpID(intenData)
        genoSnpID <- getSnpID(genoData)
        if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
        intenScanID <- getScanID(intenData)
        genoScanID <- getScanID(genoData)
        if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
  
	# what to analyze
	nSnp <- length(intenSnpID) # number of rows to analyze
        nblock <- floor(nSnp/block.size)
        remain <- nSnp-block.size*nblock
        k <- 1	# first row of first block.size
        if(remain==0) {N <- nblock} else { N <- nblock+1 }  # N is the number of blocks

	# get logical vector for samples to be excluded
        incl <- !is.element(intenScanID, scan.exclude)	
		
	# vector to store results
        mnqual <- rep(NA, nSnp)
        medqual <- rep(NA, nSnp)
	
	# read each block and calculate the mean quality over non-missing snps
	for(i in 1:N){
		if(i<=nblock) {n <- block.size}  
		else {n <- remain}
		if (verbose)
			message(paste("block number=", i, "vector index=", k))
                
		# get the data and exclude certain samples
                geno <- getGenotype(genoData, snp=c(k,n), scan=c(1,-1))
                qual <- getQuality(intenData, snp=c(k,n), scan=c(1,-1))

                # remove excluded scans
                if (length(dim(geno)) < 2) {
                  geno <- matrix(geno[incl], nrow=1)
                  qual <- matrix(qual[incl], nrow=1)
                } else {
                  geno <- geno[,incl]
                  qual <- qual[,incl]
                }
                  
                qual[is.na(geno)] <- NA
                mnqual[k:(k+n-1)] <- rowMeans(qual, na.rm=TRUE)
                medqual[k:(k+n-1)] <- rowMedians(qual, na.rm=TRUE)
		
                k <- k+n
                rm(geno); rm(qual)
	}
	qual.res <- cbind(mnqual, medqual)
        rownames(qual.res) <- intenSnpID
	colnames(qual.res) <- c("mean.quality", "median.quality")
	return(qual.res)
}

