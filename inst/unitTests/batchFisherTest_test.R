
.expFisher <- function(geno, snpID, batch, male=FALSE, CI=FALSE) {
  # expected values
  exp.batch <- sort(unique(batch))
  nbatch <- length(exp.batch)
  nsnp <- length(snpID)
  nA <- matrix(0, nrow=nsnp, ncol=nbatch)
  nB <- matrix(0, nrow=nsnp, ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- batch == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(ifelse(male, 1, 2) - geno[,thisbatch], na.rm=TRUE)
  }

  # get ci
  if (CI) {
    exp.ci1 <- rep(NA, nsnp)
    exp.ci2 <- rep(NA, nsnp)
    for (i in 1:nsnp) {
      tbl <- matrix(c(nA[i,1], nB[i,1], nA[i,2], nB[i,2]), nrow=2, ncol=2)
      try({
        exp.res <- fisher.test(tbl, conf.int=TRUE)
        ci <- exp.res[["conf.int"]]
        exp.ci1[i] <- ci[1]
        exp.ci2[i] <- ci[2]
      })
    }
    exp.ci1[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
    exp.ci2[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
    return(list(ci1=exp.ci1, ci2=exp.ci2))
  }
  
  # get pval
  if (nbatch == 2) {
    exp.pval <- rep(NA, nsnp)
    exp.or <- rep(NA, nsnp)
    for (i in 1:nsnp) {
      tbl <- matrix(c(nA[i,1], nB[i,1], nA[i,2], nB[i,2]), nrow=2, ncol=2)
      try({
        exp.res <- fisher.test(tbl, conf.int=FALSE)    
        exp.pval[i] <- exp.res[["p.value"]]   
        exp.or[i] <- exp.res[["estimate"]]
      })
    }
    exp.pval[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
    exp.or[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
    exp.or.rs <- pmin(exp.or, 1/exp.or)
    checkTrue(all(0 <= exp.or.rs[!is.na(exp.or.rs)] & exp.or.rs[!is.na(exp.or.rs)] <= 1))
    exp.ave <- rep(1/mean(exp.or.rs, na.rm=TRUE), 2)
    exp.lam <- rep(median(-2*log(exp.pval), na.rm=TRUE) / qchisq(0.5, 2), 2)
  } else {
    exp.pval <- matrix(NA, nrow=nsnp, ncol=nbatch, dimnames=list(snpID, exp.batch))
    exp.or <- matrix(NA, nrow=nsnp, ncol=nbatch, dimnames=list(snpID, exp.batch))
    for (i in 1:nsnp) {
      for (j in 1:nbatch) {
        tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                        nrow=2, ncol=2)
        try({
          exp.res <- fisher.test(tbl, conf.int=FALSE)    
          exp.pval[i,j] <- exp.res[["p.value"]]   
          exp.or[i,j] <- exp.res[["estimate"]]
        })
      }
    }
    exp.pval[pmin(rowSums(nA), rowSums(nB)) == 0,] <- NA
    exp.or[pmin(rowSums(nA), rowSums(nB)) == 0,] <- NA
    exp.ave <- rep(NA, nbatch)
    exp.lam <- rep(NA, nbatch)
    for (i in 1:nbatch) {
      exp.ave[i] <- 1/mean(pmin(exp.or[,i], 1/exp.or[,i]), na.rm=TRUE)
      exp.lam[i] <-  median(-2*log(exp.pval[,i]), na.rm=TRUE) / qchisq(0.5, 2)
    }
  }
  names(exp.ave) <- exp.batch
  names(exp.lam) <- exp.batch
  list(pval=exp.pval, or=exp.or, ave=exp.ave, lam=exp.lam)
}


check2batches <- function(nc) {
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  snpID <- getSnpID(data)
  geno <- getGenotype(data)
  exp <- .expFisher(geno, snpID, batch)

  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=FALSE)
  checkEquals(exp$ave, res$mean.or)
  checkEquals(exp$lam, res$lambda)
  checkEquals(exp$pval, res$pval[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$pval, res$pval[,2], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$or, res$oddsratio[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$or, res$oddsratio[,2], tolerance=1e-7, checkNames=FALSE)
}


check4batches <- function(nc) {
  # define batches - 4
  scanID <- 1:20
  batch <- rep(paste("batch", 1:4, sep=""), 5)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  snpID <- getSnpID(data)
  geno <- getGenotype(data)
  exp <- .expFisher(geno, snpID, batch)

  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=FALSE)
  checkEquals(exp$ave, res$mean.or)
  checkEquals(exp$lam, res$lambda)
  checkEquals(exp$pval, res$pval, tolerance=1e-7)
  checkEquals(exp$or, res$oddsratio, tolerance=1e-7)
}


checkConfInt <- function(nc) {
  # define batches - 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  snpID <- getSnpID(data)
  geno <- getGenotype(data)
  exp <- .expFisher(geno, snpID, batch, CI=TRUE)

  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=TRUE)
  checkEquals(exp$ci1, res$confint.low[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$ci2, res$confint.high[,1], tolerance=1e-7, checkNames=FALSE)
}


# test autosomes only
checkAuto <- function(nc) {
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  chr <- getChromosome(data)
  geno <- getGenotype(data)[chr <= 22,]
  snpID <- getSnpID(data)[chr <= 22]
  exp <- .expFisher(geno, snpID, batch)

  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=FALSE)
  checkEquals(exp$ave, res$mean.or)
  checkEquals(exp$lam, res$lambda)
  checkEquals(exp$pval, res$pval[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$or, res$oddsratio[,1], tolerance=1e-7, checkNames=FALSE)
}


# test X and Y chroms for males
checkXYM <- function(nc) { 
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  chr <- getChromosome(data)
  geno <- getGenotype(data)[chr %in% c(23,25), sex == "M"]
  snpID <- getSnpID(data)[chr %in% c(23,25)]
  # recode
  geno[geno == 1] <- NA
  geno[geno == 2] <- 1
  exp <- .expFisher(geno, snpID, batch[sex == "M"], male=TRUE)
  
  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=FALSE,
                         chrom.include=c(23,25), sex.include="M")
  checkEquals(exp$ave, res$mean.or)
  checkEquals(exp$lam, res$lambda)
  checkEquals(exp$pval, res$pval[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$or, res$oddsratio[,1], tolerance=1e-7, checkNames=FALSE)
}


# test X chrom for females
checkXF <- function(nc) { 
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  chr <- getChromosome(data)
  geno <- getGenotype(data)[chr == 23, sex == "F"]
  snpID <- getSnpID(data)[chr == 23]
  exp <- .expFisher(geno, snpID, batch[sex == "F"])
  
  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=FALSE,
                         chrom.include=23, sex.include="F")
  checkEquals(exp$ave, res$mean.or)
  checkEquals(exp$lam, res$lambda)
  checkEquals(exp$pval, res$pval[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$or, res$oddsratio[,1], tolerance=1e-7, checkNames=FALSE)
}


# test X and Y for males and females combined (should be error)
checkXYMF <- function(nc) {
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)
  checkException(batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, chrom.include=23))
}


checkExclude <- function(nc) {
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)
  scan.exclude <- c(1,2,10)

  # expected values
  snpID <- getSnpID(data)
  geno <- getGenotype(data)[, !(scanID %in% scan.exclude)]
  batch <- batch[!(scanID %in% scan.exclude)]
  exp <- .expFisher(geno, snpID, batch)

  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=FALSE,
                         scan.exclude=scan.exclude)
  checkEquals(exp$ave, res$mean.or)
  checkEquals(exp$lam, res$lambda)
  checkEquals(exp$pval, res$pval[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$or, res$oddsratio[,1], tolerance=1e-7, checkNames=FALSE)
}


checkSnpInclude <- function(nc) {
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)
  snpID <- getSnpID(data)
  snp.include <- sample(snpID, 10)

  # expected values
  geno <- getGenotype(data)[snpID %in% snp.include,]
  exp <- .expFisher(geno, sort(snp.include), batch)

  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=TRUE, conf.int=FALSE,
                         snp.include=snp.include)
  checkEquals(exp$ave, res$mean.or)
  checkEquals(exp$lam, res$lambda)
  checkEquals(exp$pval, res$pval[,1], tolerance=1e-7, checkNames=FALSE)
  checkEquals(exp$or, res$oddsratio[,1], tolerance=1e-7, checkNames=FALSE)
}


checkNoSnp <- function(nc) {
  # define batches - only 2
  scanID <- 1:20
  batch <- rep(paste("batch", 1:2, sep=""), 10)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # test function
  res <- batchFisherTest(data, batchVar="batch", return.by.snp=FALSE)
  checkEquals(2, length(res))
  checkIdentical(c("mean.or", "lambda"), names(res))
}


test_batchFisherTest <- function() {
  # simulate data - autosomes only
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=22,
                           n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)

  # run tests
  check2batches(nc)
  check4batches(nc)
  checkConfInt(nc)
  checkExclude(nc)
  checkSnpInclude(nc)
  checkNoSnp(nc)

  # clean up
  close(nc)
  file.remove(ncfile)

  
  # simulate data - 26 chroms
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                           n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)

  # run tests
  checkAuto(nc)
  checkXYM(nc)
  checkXF(nc)
  checkXYMF(nc)

  # clean up
  close(nc)
  file.remove(ncfile)
}
