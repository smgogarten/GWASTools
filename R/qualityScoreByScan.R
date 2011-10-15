qualityScoreByScan <- 
function(intenData,
         genoData,
         snp.exclude=NULL,
         verbose=TRUE) 
{

	# get mean and median quality score over snps within a sample for non-missing genotype calls
	# exclude Y chromosome SNPs for females

        # check that intenData has quality
        if (!hasQuality(intenData)) stop("quality not found in intenData")
  
        # check that dimensions of intenData and genoData are equal
        intenSnpID <- getSnpID(intenData)
        genoSnpID <- getSnpID(genoData)
        if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
        intenScanID <- getScanID(intenData)
        genoScanID <- getScanID(genoData)
        if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")

        # check that sex is present in annotation
        if (hasSex(intenData)) {
          sex <- getSex(intenData)
        } else if (hasSex(genoData)) {
          sex <- getSex(genoData)
        } else stop("sex not found in intenData or genoData")

 	fem <- which(sex == "F")

	# define a logical indicator for Y chromosome
	chr <- getChromosome(intenData, char=TRUE)
	ychr <- chr == "Y"

        # logical vector for SNPs to keep
        incl <- !is.element(intenSnpID, snp.exclude)

        #define empty vectors to fill in with mean and median genotype confidences
	nScan <- length(intenScanID)
	mnqual <- rep(NA, nScan); names(mnqual) <- intenScanID
	medqual <- rep(NA, nScan); names(medqual) <- intenScanID

        # for each sample
	for(i in 1:nScan)
	{
		if (verbose & (i %% 100 == 0))
			message(paste("scan", i, "of", nScan))

                geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))
                qual <- getQuality(intenData, snp=c(1,-1), scan=c(i,1))
		# omit Y chrom from genotypes and qualities if the sample is from a female
		if(is.element(i,fem)){
                        geno <- geno[!ychr & incl]
                        qual <- qual[!ychr & incl]
                } else {
                        geno <- geno[incl]
                        qual <- qual[incl]
                }
                
		# take the mean and median genotype quality score from the nonmissing genotypes in the sample	
                mnqual[i] <- mean(qual[!is.na(geno)], na.rm=TRUE)
                medqual[i] <- median(qual[!is.na(geno)], na.rm=TRUE)
                rm(geno); rm(qual)
	}       
	qual.res <- cbind(mnqual, medqual)
        rownames(qual.res) <- intenScanID
	colnames(qual.res) <- c("mean.quality", "median.quality")
	return(qual.res)
}
