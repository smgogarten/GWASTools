
#############################################################################
#
# Batch Chi-square calculations
#
#  genoData -- object of type GenotypeData
#  batchVar	-- sample batch (name of Annotation column)
#  correct      -- apply Yates continuity correction
#  verbose	-- whether to show the progress
#
# Return: average chi-square values for batches
#
#############################################################################

.batchSnpInclude <- function(genoData, snp.include, chrom.include) {
    if (!is.null(snp.include)) {
        getSnpID(genoData) %in% snp.include
    } else {
        getChromosome(genoData) %in% chrom.include
    }
}

.batchSexCheck <- function(genoData, snpIndex, sex.include) {
    stopifnot(hasSex(genoData))
    chrom.include <- unique(getChromosome(genoData, index=snpIndex))
    if (any(c(XchromCode(genoData), YchromCode(genoData)) %in% chrom.include) &
        length(sex.include) > 1) {
        stop("sex chromosome batch tests cannot be done for both sexes combined")
    }
}

.getBatch <- function(genoData, batchVar, scan.exclude) {
    stopifnot(is.character(batchVar))
    stopifnot(hasScanVariable(genoData, batchVar))
    
    batch <- getScanVariable(genoData, batchVar)

    ## exclude scans
    if (!is.null(scan.exclude)) {
        batch[getScanID(genoData) %in% scan.exclude] <- NA
    }
    
    id <- levels(factor(batch))
    nBatch <- length(id)
    if (nBatch < 2) stop(batchVar, " must have at least two unique values")
    list(var=batch, id=id, n=nBatch)
}

.batchGenotypeCounts <- function(genoData, snpIndex, batch, sex.include, verbose) {
    if (verbose)
        message(date(), "\tLoad genotype dataset ...")
    
    chrom <- getChromosome(genoData)
    xchr <- chrom[snpIndex] == XchromCode(genoData)
    ychr <- chrom[snpIndex] == YchromCode(genoData)
    snpID <- getSnpID(genoData, index=snpIndex)
    nSnp <- length(snpID)
    
    ## the number of A and B alleles for plates and SNPs
    nA <- matrix(0, nrow=nSnp, ncol=batch$n)
    nB <- matrix(0, nrow=nSnp, ncol=batch$n)
    
    sex <- getSex(genoData)
    nSamp <- nscan(genoData)
    for (iSamp in 1:nSamp) {
        if (!is.na(batch$var[iSamp]) & (sex[iSamp] %in% sex.include)) {
            i <- match(batch$var[iSamp], batch$id)
            geno <- getGenotype(genoData, snp=c(1,-1), scan=c(iSamp,1))
            geno <- geno[snpIndex]
            ## if subject is male, recode X & Y genotypes from (0,1,2) to (0,NA,1)
            ## since we are treating missing genotypes as 0 for counting purposes,
            ## code this as (0,0,1)
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
    ## allele counts per snp
    nAlleles <- matrix(c(rowSums(nA), rowSums(nB)), nrow=nSnp, ncol=2,
                       dimnames=list(snpID, c("nA","nB")))
    ## minor allele counts
    minorCounts <- rowMin(nAlleles)

    list(nA=nA, nB=nB, nAlleles=nAlleles, minor=minorCounts)
}

## minimum expected frequency in each cell
.minExpFreq <- function(n1A, n1B, n2A, n2B) {

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
                           snp.include = NULL,
                           chrom.include = 1:22,
                           sex.include = c("M", "F"),
                           scan.exclude = NULL,
                           return.by.snp = FALSE,
                           correct = TRUE,
                           verbose = TRUE) {

    snpIndex <- .batchSnpInclude(genoData, snp.include, chrom.include)
        
    .batchSexCheck(genoData, snpIndex, sex.include)

    batch <- .getBatch(genoData, batchVar, scan.exclude) 

    genoCounts <- .batchGenotypeCounts(genoData, snpIndex, batch, sex.include, verbose)

    ## prepare data
    ave <- setNames(double(batch$n), batch$id)
    lambda <- setNames(double(batch$n), batch$id)
    if (return.by.snp) {
        snpID <- rownames(genoCounts$nAlleles)
        nSnp <- length(snpID)
        Chi.Batch <- matrix(nrow=nSnp, ncol=batch$n, dimnames=list(snpID, batch$id))
        Exp.Batch <- matrix(nrow=nSnp, ncol=batch$n, dimnames=list(snpID, batch$id))
    }

    ## chi-square, vector calculation
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

        ## don't forget "1.0" to avoid integer overflow
        ntot <- (n1A + n1B + n2A + n2B)
        if (correct) {
            yates <- pmax(0, abs(n1A*n2B - n1B*n2A) - 0.5*ntot)
            rv <- 1.0 * ntot * (yates)^2 /
                (1.0 * (n1A+n1B) * (n1A+n2A) * (n2B+n1B) * (n2B+n2A) )
        } else {
            rv <- 1.0 * ntot * (n1A*n2B - n1B*n2A)^2 /
                (1.0 * (n1A+n1B) * (n1A+n2A) * (n2B+n1B) * (n2B+n2A) )
        }
        rv[!is.finite(rv)] <- NA

        ## return NA for intensity-only and monomorphic SNPs
        rv[genoCounts$minor == 0] <- NA

        ## save
        if (return.by.snp) {
            Chi.Batch[,i] <- rv
            Exp.Batch[,i] <- .minExpFreq(n1A, n1B, n2A, n2B)
        }

        ## mean chisq
        ave[i] <- mean(rv, na.rm=TRUE)
        ## genomic inflation factor
        lambda[i] <- median(rv, na.rm=TRUE) / qchisq(0.5, 1)

        if (verbose)
            message(date(), "\t", i, "/", batch$n, "\t", batch$id[i])

        ## if there are only 2 batches, no need to repeat calculation twice
        if (batch$n == 2) {
            ave[i+1] <- ave[i]
            lambda[i+1] <- lambda[i]
            if (return.by.snp) {
                Chi.Batch[,i+1] <- Chi.Batch[,i]
                Exp.Batch[,i+1] <- Exp.Batch[,i]
            }
            break
        }
    }

    if (return.by.snp) {
        list("mean.chisq"=ave, "lambda"=lambda, "chisq"=Chi.Batch, "allele.counts"=genoCounts$nAlleles, "min.exp.freq"=Exp.Batch)
    } else {
        list("mean.chisq"=ave, "lambda"=lambda)
    }
}

