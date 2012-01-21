
#############################################################################
#
# Batch Fisher Exact Test calculations
#
#  genoData -- object of type GenotypeData
#  batchVar	-- sample batch (name of Annotation column)
#  verbose	-- whether to show the progress
#
# Return: average p values for batches
#
#############################################################################
batchFisherTest <- function(genoData,
                             batchVar, 
                             chrom.include = 1:22,
                             sex.include = c("M", "F"),
                           scan.exclude = NULL,
                           return.by.snp = FALSE,
                             conf.int = FALSE,
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

        chrom <- getChromosome(genoData)
        chromIndex <- chrom %in% chrom.include
        xchr <- chrom[chromIndex] == XchromCode(genoData)
        ychr <- chrom[chromIndex] == YchromCode(genoData)
        snpID <- getSnpID(genoData, index=chromIndex)
        nSnp <- length(snpID)
        
	# prepare data
        ave <- double(nBatch)
        names(ave) <- id
        lambda <- double(nBatch)
        names(lambda) <- id
        if (return.by.snp) {
          Pval.Batch <- matrix(nrow=nSnp, ncol=nBatch, dimnames=list(snpID, id))
          OR.Batch <- matrix(nrow=nSnp, ncol=nBatch, dimnames=list(snpID, id))
          Exp.Batch <- matrix(nrow=nSnp, ncol=nBatch, dimnames=list(snpID, id))
          if (conf.int) {
            CI1.Batch <- matrix(nrow=nSnp, ncol=nBatch, dimnames=list(snpID, id))
            CI2.Batch <- matrix(nrow=nSnp, ncol=nBatch, dimnames=list(snpID, id))
          }
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

	for (i in 1:nBatch)
	{
		n1A <- nA[, i]; n1B <- nB[, i]
		if (length(id) > 2)
		{
			n2A <- rowSums(nA[, -i]); n2B <- rowSums(nB[, -i])
		} else {
			n2A <- nA[, -i]; n2B <- nB[, -i]
		}

                # loop over SNPs
                pval <- rep(NA, nSnp)
                or <- rep(NA, nSnp)
                if (conf.int) {
                  ci <- matrix(nrow=nSnp, ncol=2)
                }
                for (j in 1:nSnp) {
                  # skip intensity-only and monomorphic SNPs (return NA)
                  if (minorCounts[j] > 0) {
                    tbl <- matrix(c(n1A[j], n1B[j], n2A[j], n2B[j]), nrow=2, ncol=2)
                    if (all(!is.na(tbl))) {
                      # try fisher's exact test, but don't quit on error
                      # result will be NA
                      try({
                        rv <- fisher.test(tbl, conf.int=conf.int)
                        pval[j] <- rv$p.value
                        or[j] <- rv$estimate
                        if (conf.int) {
                          ci[j,] <- rv$conf.int
                        }
                      })
                    }
                  }
                }
                
		# save
                if (return.by.snp) {
                  Pval.Batch[,i] <- pval
                  OR.Batch[,i] <- or
                  if (conf.int) {
                    CI1.Batch[,i] <- ci[,1]
                    CI2.Batch[,i] <- ci[,2]
                  }
                  Exp.Batch[,i] <- GWASTools:::minExpFreq(n1A, n1B, n2A, n2B)
                }
                
                # odds ratio goes from 0 to infinity
                # rescale so it is between 0 and 1 for computing mean
                # then take inverse so mean odds ratios are > 1
                ave[i] <- 1/mean(pmin(or, 1/or), na.rm=TRUE)
  	        # genomic inflation factor
		lambda[i] <- median(-2*log(pval), na.rm=TRUE) / 1.39

		if (verbose)
			message(date(), "\t", i, "/", nBatch, "\t", id[i])

                # if there are only 2 batches, no need to repeat calculation twice
                if (nBatch == 2) {               
                  ave[i+1] <- ave[i]
                  lambda[i+1] <- lambda[i]
                  if (return.by.snp) {
                    Pval.Batch[,i+1] <- Pval.Batch[,i]
                    OR.Batch[,i+1] <- OR.Batch[,i]
                    if (conf.int) {
                      CI1.Batch[,i+1] <- CI1.Batch[,i]
                      CI2.Batch[,i+1] <- CI2.Batch[,i]
                    }
                    Exp.Batch[,i+1] <- Exp.Batch[,i]
                  }
                  break
                }
	}

        if (return.by.snp) {
          if (conf.int) {
            fisher.batch <- list("mean.or"=ave, "lambda"=lambda, "pval"=Pval.Batch, "oddsratio"=OR.Batch, "confint.low"=CI1.Batch, "confint.high"=CI2.Batch, "allele.counts"=nAlleles, "min.exp.freq"=Exp.Batch)
          } else {
            fisher.batch <- list("mean.or"=ave, "lambda"=lambda, "pval"=Pval.Batch, "oddsratio"=OR.Batch, "allele.counts"=nAlleles, "min.exp.freq"=Exp.Batch)
          }
        } else {
          fisher.batch <- list("mean.or"=ave, "lambda"=lambda)
        }
        
	if (verbose)
		message(date(), "\t", "End ...")

        if (!is.null(outfile)) {
          fileOut <- paste(outfile, "RData", sep=".")
          save(fisher.batch, file=fileOut, compress=TRUE)

          # save the warnings
          warn <- warnings()
          if (!is.null(warn)) {
            warnfileOut <- paste(outfile, "warnings", "RData", sep=".")
            save(warn, file=warnfileOut)
          }
          return(invisible(NULL))
        } else {
          return(fisher.batch)
        }
}

