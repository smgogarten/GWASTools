
#############################################################################
#
# Batch Chi-square calculations
#
#  genoData -- object of type GenotypeData
#  batchVar	-- sample batch (name of Annotation column)
#  correct      -- apply Yates continuity correction
#  verbose	-- whether to show the progress
#  outfile	-- name of a file to save output to
#
# Return: average chi-square values for batches
#
#############################################################################

# minimum expected frequency in each cell
minExpFreq <- function(n1A, n1B, n2A, n2B) {

        ntot <- (n1A + n1B + n2A + n2B)
        f1 <- n1A + n1B
        f2 <- n2A + n2B
        fA <- n1A + n2A
        fB <- n1B + n2B
        exp1A <- (f1 * fA) / ntot
        exp2A <- (f2 * fA) / ntot
        exp1B <- (f1 * fB) / ntot
        exp2B <- (f2 * fB) / ntot
        minExp <- pmin(exp1A, exp2A, exp1B, exp2B)
}

batchChisqTest <- function(genoData,
                             batchVar, 
                             chrom.include = 1:22,
                             sex.include = c("M", "F"),
                           scan.exclude = NULL,
                           return.by.snp = FALSE,
                             correct = TRUE,
                             verbose = TRUE,
                           outfile = NULL)
{


	# check
	stopifnot(is.character(batchVar))
        stopifnot(hasScanVariable(genoData, batchVar))
        stopifnot(hasSex(genoData))
	stopifnot(is.null(outfile) | is.character(outfile))

        if (any(c(XchromCode(genoData), YchromCode(genoData)) %in% chrom.include) &
            length(sex.include) > 1) {
          stop("sex chromosome batch tests cannot be done for both sexes combined")
        }
        
        batch <- getScanVariable(genoData, batchVar)
        sex <- getSex(genoData)
        scanID <- getScanID(genoData)

        # exclude scans
        if (!is.null(scan.exclude)) {
          batch[scanID %in% scan.exclude] <- NA
        }
        
	id <- levels(factor(batch))
        nBatch <- length(id)
	if (nBatch < 2)
		stop("The level of batch should be >= 2!")

         # if nBatch is 2, only compute and store results once
         if (nBatch == 2) {
           outBatch <- 1
         } else {
           outBatch <- nBatch
         }
                
        chrom <- getChromosome(genoData)
        chromIndex <- chrom %in% chrom.include
        xchr <- chrom[chromIndex] == XchromCode(genoData)
        ychr <- chrom[chromIndex] == YchromCode(genoData)
        snpID <- getSnpID(genoData, index=chromIndex)
        nSnp <- length(snpID)
        
	# prepare data
        ave <- double(outBatch)
        names(ave) <- id[1:outBatch]
        lambda <- double(outBatch)
        names(lambda) <- id[1:outBatch]
        if (return.by.snp) {
          Chi.Batch <- matrix(nrow=nSnp, ncol=outBatch, dimnames=list(snpID, id[1:outBatch]))
          Exp.Batch <- matrix(nrow=nSnp, ncol=outBatch, dimnames=list(snpID, id[1:outBatch]))
        }

	# the number of A alleles for plates and SNPs
	nA <- matrix(0, nrow=nSnp, ncol=nBatch)
	# the number of B alleles for plates and SNPs
	nB <- matrix(0, nrow=nSnp, ncol=nBatch)

	if (verbose)
		message(date(), "\tLoad genotype dataset ...")
        nSamp <- nscan(genoData)
	for (iSamp in 1:nSamp)
	{
		if (!is.na(batch[iSamp]) & (sex[iSamp] %in% sex.include))
		{
			i <- match(batch[iSamp], id)
                        geno <- getGenotype(genoData, snp=c(1,-1), scan=c(iSamp,1))
                        geno <- geno[chromIndex]
                        # if subject is male, recode X & Y genotypes from (0,1,2) to (0,NA,1)
                        # since we are treating missing genotypes as 0 for counting purposes,
                        # code this as (0,0,1)
			gA <- geno
                        gA[is.na(geno)] <- 0
                        if (sex[iSamp] == "M") {
                          gA[(xchr | ychr) & gA == 1] <- 0
                          gA[(xchr | ychr) & gA == 2] <- 1
                        }
			nA[, i] <- nA[, i] + gA
                        
			gB <- 2-geno
                        gB[is.na(geno)] <- 0
                        if (sex[iSamp] == "M") {
                          gB[(xchr | ychr) & gB == 1] <- 0
                          gB[(xchr | ychr) & gB == 2] <- 1
                        }
			nB[, i] <- nB[, i] + gB

			if (verbose & (iSamp %% 500 == 0))
				message(date(), "\t\t", iSamp, "/", nSamp)
		}
	}
        # allele counts per snp
        nAlleles <- matrix(c(rowSums(nA), rowSums(nB)), nrow=nSnp, ncol=2,
                           dimnames=list(snpID, c("nA","nB")))
        # minor allele counts
        minorCounts <- rowMin(nAlleles)
        
	if (verbose)
		message(date(), "\tStart ...")

	# chi-square, vector calculation

	for (i in 1:nBatch)
	{
		n1A <- nA[, i]; n1B <- nB[, i]
		if (length(id) > 2)
		{
			n2A <- rowSums(nA[, -i]); n2B <- rowSums(nB[, -i])
		} else {
			n2A <- nA[, -i]; n2B <- nB[, -i]
		}

		# don't forget "1.0" to avoid integer overflow
                ntot <- (n1A + n1B + n2A + n2B)
                if (correct) {
                  rv <- 1.0 * ntot * (abs(n1A*n2B - n1B*n2A) - 0.5*ntot)^2 /
                    (1.0 * (n1A+n1B) * (n1A+n2A) * (n2B+n1B) * (n2B+n2A) )
                } else {
                  rv <- 1.0 * ntot * (n1A*n2B - n1B*n2A)^2 /
                    (1.0 * (n1A+n1B) * (n1A+n2A) * (n2B+n1B) * (n2B+n2A) )
                }
		rv[!is.finite(rv)] <- NA

                # return NA for intensity-only and monomorphic SNPs
                rv[minorCounts == 0] <- NA

		# save
                if (return.by.snp) {
                  Chi.Batch[,i] <- rv
                  Exp.Batch[,i] <- GWASTools:::minExpFreq(n1A, n1B, n2A, n2B)
                }

  	        # genomic inflation factor
                ave[i] <- mean(rv, na.rm=TRUE)
                lambda[i] <- median(rv, na.rm=TRUE) / 0.456

		if (verbose)
			message(date(), "\t", i, "/", nBatch, "\t", id[i])
                
                # if there are only 2 batches, no need to repeat calculation twice
                if (nBatch == 2) {
                  break
                }
	}
        
        if (return.by.snp) {
          res <- list("mean.chisq"=ave, "lambda"=lambda, "chisq"=Chi.Batch, "allele.counts"=nAlleles, "min.exp.freq"=Exp.Batch)
        } else {
          res <- list("mean.chisq"=ave, "lambda"=lambda)
        }

	if (verbose)
		message(date(), "\t", "End ...")

        if (!is.null(outfile)) {
          fileOut <- paste(outfile, "RData", sep=".")
          save(res, file=fileOut, compress=TRUE)

          # save the warnings
          warn <- warnings()
          if (!is.null(warn)) {
            warnfileOut <- paste(outfile, "warnings", "RData", sep=".")
            save(warn, file=warnfileOut)
          }
          return(invisible(NULL))
        } else {
          return(res)
        }

}

