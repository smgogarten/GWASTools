missingGenotypeByScanChrom <- 
function(
		genoData, # object of type GenotypeData
		snp.exclude=NULL,
		verbose=TRUE) 
{

	# get start and count row indices for each chromosome type
        chrom <- getChromosome(genoData, char=TRUE)
        uniqChrom <- unique(chrom)
        nChrom <- length(uniqChrom)
	indices <- matrix(NA, nChrom, 2, dimnames=list(uniqChrom, c("start","stop")))
	for(i in 1:nChrom) indices[i,] <- range(which(is.element(chrom, uniqChrom[i])))
	# get snps per chromosome
	spc <- apply(indices, 1, function(x) x[2]-x[1]+1)

	# logical vectors for excluding snps by chromosome
        snpID <- getSnpID(genoData)
	if(length(snp.exclude)!=0) {
		exclude <- !is.element(snpID, snp.exclude)
		chrex <- vector("list", nChrom)
		for(j in 1:nChrom){ 
			chrex[[j]] <- exclude[indices[j,"start"]:indices[j,"stop"]]
			spc[j] <- sum(chrex[[j]]) # number of non-excluded snps per chromosome
		}
	}

        scanID <- getScanID(genoData)
        nScan <- length(scanID)

	# matrix to hold missing counts
	miss.cnt <- matrix(NA, length(scanID), nChrom, dimnames=list(scanID, uniqChrom)) 

	for(i in 1:nScan)
	{	# for each sample
		if (verbose & (i %% 100 == 0))
			message(paste("scan",i,"of",nScan))
		# get all snps for sample i
                geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))

		for(j in 1:nChrom) {	# for each chromosome
			x <- geno[indices[j,"start"]:indices[j,"stop"]]	# genotypes for the jth chromosome
			if(length(snp.exclude)!=0) x <- x[chrex[[j]]]	# remove snps to be excluded
			miss.cnt[i,j] <- length(x[is.na(x)])
			rm(x)
		}
		rm(geno)
	}

        miss.frac <- rowSums(miss.cnt) / sum(spc)
        sex <- getSex(genoData)
        # for females, exclude Y chromosome
        female <- which(sex == "F")
        notY <- uniqChrom != "Y"
        miss.frac[female] <- rowSums(miss.cnt[female, notY, drop=FALSE]) / sum(spc[notY])
        
	miss <- list(miss.cnt, spc, miss.frac)
        names(miss) <- c("missing.counts", "snps.per.chr", "missing.fraction")
	return(miss)
}

