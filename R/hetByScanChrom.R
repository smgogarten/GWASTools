


### returns heterozygosity over samples for each chromosome, plus autosomal average heterozygosity over chromosomes (NOT SNPS)


hetByScanChrom <- function(
                        genoData, # object of type GenotypeData
			snp.exclude = NULL, 
			verbose = TRUE) 
{
    # get start and count row indices for each chromosome type
    chrom <- getChromosome(genoData, char=TRUE)
    uniqChrom <- unique(chrom)
    nChrom <- length(uniqChrom)
    indices <- matrix(NA, nChrom, 2, dimnames = list(uniqChrom, 
        c("start", "stop")))
    for (i in 1:nChrom)
      indices[i, ] <- range(which(is.element(chrom, uniqChrom[i])))

    # logical vectors for excluding snps by chromosome
    snpID <- getSnpID(genoData)
    if (length(snp.exclude) != 0) {
        exclude <- !is.element(snpID, snp.exclude)
        chrex <- vector("list", nChrom)
        for (j in 1:nChrom) {
            chrex[[j]] <- exclude[indices[j, "start"]:indices[j, "stop"]]
        }
    }

    scanID <- getScanID(genoData)
    nScan <- length(scanID)
    
    # matrix to hold heterozygote counts
    het <- matrix(NA, nScan, nChrom, dimnames = list(scanID, uniqChrom))
    # matrix to hold non-missing genotype counts
    nm <- matrix(NA, nScan, nChrom, dimnames = list(scanID, uniqChrom))
    # vector to hold autosomal het
    A <- rep(NA, nScan)

    # for each sample
    for (i in 1:nScan) { 
        if (verbose & (i %% 100 == 0)) 
            message(paste("scan", i, "of", nScan))
	# get all snps for sample i
        geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))

        # for each chromosome
        for (j in 1:nChrom) {
            # genotypes for the jth chromosome
            x <- geno[indices[j, "start"]:indices[j, "stop"]]
            # remove snps to be excluded
            if(length(snp.exclude) != 0) x <- x[chrex[[j]]]
            nm[i, j] <- length(x[!is.na(x)])
            het[i, j] <- length(x[!is.na(x) & x == 1])
            rm(x)
        }
        rm(geno)
	# heterozygous fraction for autosomes
        index <- is.element(uniqChrom, as.character(autosomeCode(genoData)))
        A[i] <- sum(het[i, index])/sum(nm[i, index])
    }
    het <- het/nm  # heterozygous fraction
    het <- cbind(het, A)
    return(het)
}
