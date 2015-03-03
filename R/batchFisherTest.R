
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
                            snp.include = NULL,
                            chrom.include = 1:22,
                            sex.include = c("M", "F"),
                            scan.exclude = NULL,
                            return.by.snp = FALSE,
                            conf.int = FALSE,
                            verbose = TRUE) {

    snpIndex <- .batchSnpInclude(genoData, snp.include, chrom.include)
        
    .batchSexCheck(genoData, snpIndex, sex.include)

    batch <- .getBatch(genoData, batchVar, scan.exclude) 

    genoCounts <- .batchGenotypeCounts(genoData, snpIndex, batch, sex.include, verbose)
    nSnp <- nrow(genoCounts$nAlleles)

    ## prepare data
    ave <- setNames(double(batch$n), batch$id)
    lambda <- setNames(double(batch$n), batch$id)
    if (return.by.snp) {
        snpID <- rownames(genoCounts$nAlleles)
        Pval.Batch <- matrix(nrow=nSnp, ncol=batch$n, dimnames=list(snpID, batch$id))
        OR.Batch <- matrix(nrow=nSnp, ncol=batch$n, dimnames=list(snpID, batch$id))
        Exp.Batch <- matrix(nrow=nSnp, ncol=batch$n, dimnames=list(snpID, batch$id))
        if (conf.int) {
            CI1.Batch <- matrix(nrow=nSnp, ncol=batch$n, dimnames=list(snpID, batch$id))
            CI2.Batch <- matrix(nrow=nSnp, ncol=batch$n, dimnames=list(snpID, batch$id))
        }
    }

    for (i in 1:batch$n) {
        n1A <- genoCounts$nA[, i]
        n1B <- genoCounts$nB[, i]
        if (length(batch$id) > 2) {
            n2A <- rowSums(genoCounts$nA[, -i])
            n2B <- rowSums(genoCounts$nB[, -i])
        } else {
            n2A <- genoCounts$nA[, -i]
            n2B <- genoCounts$nB[, -i]
        }

        ## loop over SNPs
        pval <- rep(NA, nSnp)
        or <- rep(NA, nSnp)
        if (conf.int) {
            ci <- matrix(nrow=nSnp, ncol=2)
        }
        for (j in 1:nSnp) {
            ## skip intensity-only and monomorphic SNPs (return NA)
            if (genoCounts$minor[j] > 0) {
                tbl <- matrix(c(n1A[j], n1B[j], n2A[j], n2B[j]), nrow=2, ncol=2)
                if (all(!is.na(tbl))) {
                    ## try fisher's exact test, but don't quit on error
                    ## result will be NA
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

        ## save
        if (return.by.snp) {
            Pval.Batch[,i] <- pval
            OR.Batch[,i] <- or
            if (conf.int) {
                CI1.Batch[,i] <- ci[,1]
                CI2.Batch[,i] <- ci[,2]
            }
            Exp.Batch[,i] <- .minExpFreq(n1A, n1B, n2A, n2B)
        }

        ## odds ratio goes from 0 to infinity
        ## rescale so it is between 0 and 1 for computing mean
        ## then take inverse so mean odds ratios are > 1
        ave[i] <- 1/mean(pmin(or, 1/or), na.rm=TRUE)
        ## genomic inflation factor
        lambda[i] <- median(-2*log(pval), na.rm=TRUE) / qchisq(0.5, 2)

        if (verbose)
            message(date(), "\t", i, "/", batch$n, "\t", batch$id[i])

        ## if there are only 2 batches, no need to repeat calculation twice
        if (batch$n == 2) {
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
            list("mean.or"=ave, "lambda"=lambda, "pval"=Pval.Batch, "oddsratio"=OR.Batch, "confint.low"=CI1.Batch, "confint.high"=CI2.Batch, "allele.counts"=genoCounts$nAlleles, "min.exp.freq"=Exp.Batch)
        } else {
            list("mean.or"=ave, "lambda"=lambda, "pval"=Pval.Batch, "oddsratio"=OR.Batch, "allele.counts"=genoCounts$nAlleles, "min.exp.freq"=Exp.Batch)
        }
    } else {
        list("mean.or"=ave, "lambda"=lambda)
    }
}

