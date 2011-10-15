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
  exp.batch <- unique(batch)
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  nB <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- scanAnnot$batch == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(2-geno[,thisbatch], na.rm=TRUE)
  }

  # get chisq
  exp.chisq <- matrix(NA, nrow=nsnp(data), ncol=1, dimnames=list(snpID, exp.batch[1]))
  for (i in 1:nsnp(data)) {
    tbl <- matrix(c(nA[i,1], nB[i,1], nA[i,2], nB[i,2]), nrow=2, ncol=2)
    try({
      chisq <- chisq.test(tbl, correct=TRUE)
      exp.chisq[i] <- chisq[["statistic"]]
    })
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- mean(exp.chisq, na.rm=TRUE)
  names(exp.ave) <- exp.batch[1]
  exp.lam <- median(exp.chisq, na.rm=TRUE) / 0.456
  names(exp.lam) <- exp.batch[1]

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=TRUE)
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)
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
  exp.batch <- unique(batch)
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  nB <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- scanAnnot$batch == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(2-geno[,thisbatch], na.rm=TRUE)
  }
  exp.chisq <- matrix(NA, nrow=nsnp(data), ncol=nbatch, dimnames=list(snpID, exp.batch))
  for (i in 1:nsnp(data)) {
    for (j in 1:nbatch) {
      tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                      nrow=2, ncol=2)
      try({
        chisq <- chisq.test(tbl, correct=TRUE)
        exp.chisq[i,j] <- chisq[["statistic"]]
      })
    }
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- rep(NA, nbatch)
  names(exp.ave) <- exp.batch
  exp.lam <- rep(NA, nbatch)
  names(exp.lam) <- exp.batch
  for (i in 1:nbatch) {
    exp.ave[i] <- mean(exp.chisq[,i], na.rm=TRUE)
    exp.lam[i] <- median(exp.chisq[,i], na.rm=TRUE) / 0.456
  }

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=TRUE)
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)

}


# test without Yates correction
checkNoYates <- function(nc) {
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
  exp.batch <- unique(batch)
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  nB <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- scanAnnot$batch == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(2-geno[,thisbatch], na.rm=TRUE)
  }
  exp.chisq <- matrix(NA, nrow=nsnp(data), ncol=nbatch, dimnames=list(snpID, exp.batch))
  for (i in 1:nsnp(data)) {
    for (j in 1:nbatch) {
      tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                      nrow=2, ncol=2)
      try({
        chisq <- chisq.test(tbl, correct=FALSE)
        exp.chisq[i,j] <- chisq[["statistic"]]
      })
    }
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- rep(NA, nbatch)
  names(exp.ave) <- exp.batch
  exp.lam <- rep(NA, nbatch)
  names(exp.lam) <- exp.batch
  for (i in 1:nbatch) {
    exp.ave[i] <- mean(exp.chisq[,i], na.rm=TRUE)
    exp.lam[i] <- median(exp.chisq[,i], na.rm=TRUE) / 0.456
  }

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=TRUE, correct=FALSE)
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)
}


# test autosomes only
checkAuto <- function(nc) {
  # define batches - 4
  scanID <- 1:20
  batch <- rep(paste("batch", 1:4, sep=""), 5)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  chr <- getChromosome(data)
  geno <- getGenotype(data)[chr <= 22,]
  snpID <- getSnpID(data)[chr <= 22]
  exp.batch <- unique(batch)
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow(geno), ncol=nbatch)
  nB <- matrix(0, nrow(geno), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- scanAnnot$batch == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(2-geno[,thisbatch], na.rm=TRUE)
  }

  # get chisq
  exp.chisq <- matrix(NA, nrow=nrow(geno), ncol=nbatch, dimnames=list(snpID, exp.batch))
  for (i in 1:nrow(geno)) {
    for (j in 1:nbatch) {
      tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                      nrow=2, ncol=2)
      try({
        chisq <- chisq.test(tbl, correct=TRUE)
        exp.chisq[i,j] <- chisq[["statistic"]]
      })
    }
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- rep(NA, nbatch)
  names(exp.ave) <- exp.batch
  exp.lam <- rep(NA, nbatch)
  names(exp.lam) <- exp.batch
  for (i in 1:nbatch) {
    exp.ave[i] <- mean(exp.chisq[,i], na.rm=TRUE)
    exp.lam[i] <- median(exp.chisq[,i], na.rm=TRUE) / 0.456
  }

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=TRUE)
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)
}


# test X and Y chroms for males
checkXYM <- function(nc) { 
  # define batches - 4
  scanID <- 1:20
  batch <- rep(paste("batch", 1:4, sep=""), 5)
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
  
  exp.batch <- unique(batch)
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow(geno), ncol=nbatch)
  nB <- matrix(0, nrow(geno), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- scanAnnot$batch[sex == "M"] == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(1-geno[,thisbatch], na.rm=TRUE)
  }

  # get chisq
  exp.chisq <- matrix(NA, nrow=nrow(geno), ncol=nbatch, dimnames=list(snpID, exp.batch))
  for (i in 1:nrow(geno)) {
    for (j in 1:nbatch) {
      tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                      nrow=2, ncol=2)
      try({
        chisq <- chisq.test(tbl, correct=TRUE)
        exp.chisq[i,j] <- chisq[["statistic"]]
      })
    }
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- rep(NA, nbatch)
  names(exp.ave) <- exp.batch
  exp.lam <- rep(NA, nbatch)
  names(exp.lam) <- exp.batch
  for (i in 1:nbatch) {
    exp.ave[i] <- mean(exp.chisq[,i], na.rm=TRUE)
    exp.lam[i] <- median(exp.chisq[,i], na.rm=TRUE) / 0.456
  }

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=TRUE,
                        chrom.include=c(23,25), sex.include="M")
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)
}


# test X chrom for females
checkXF <- function(nc) { 
  # define batches - 4
  scanID <- 1:20
  batch <- rep(paste("batch", 1:4, sep=""), 5)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # expected values
  chr <- getChromosome(data)
  geno <- getGenotype(data)[chr == 23, sex == "F"]
  snpID <- getSnpID(data)[chr == 23]
  
  exp.batch <- unique(batch)
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow(geno), ncol=nbatch)
  nB <- matrix(0, nrow(geno), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- scanAnnot$batch[sex == "F"] == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(2-geno[,thisbatch], na.rm=TRUE)
  }

  # get chisq
  exp.chisq <- matrix(NA, nrow=nrow(geno), ncol=nbatch, dimnames=list(snpID, exp.batch))
  for (i in 1:nrow(geno)) {
    for (j in 1:nbatch) {
      tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                      nrow=2, ncol=2)
      try({
        chisq <- chisq.test(tbl, correct=TRUE)
        exp.chisq[i,j] <- chisq[["statistic"]]
      })
    }
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- rep(NA, nbatch)
  names(exp.ave) <- exp.batch
  exp.lam <- rep(NA, nbatch)
  names(exp.lam) <- exp.batch
  for (i in 1:nbatch) {
    exp.ave[i] <- mean(exp.chisq[,i], na.rm=TRUE)
    exp.lam[i] <- median(exp.chisq[,i], na.rm=TRUE) / 0.456
  }

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=TRUE,
                        chrom.include=23, sex.include="F")
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)
}


# test X and Y for males and females combined (should be error)
checkXYMF <- function(nc){
  # define batches - 4
  scanID <- 1:20
  batch <- rep(paste("batch", 1:4, sep=""), 5)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)
  checkException(batchChisqTest(data, batchVar="batch", chrom.include=23))
}


checkFileOut <- function(nc) {
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
  exp.batch <- unique(batch)
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  nB <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- scanAnnot$batch == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(2-geno[,thisbatch], na.rm=TRUE)
  }

  # get chisq
  exp.chisq <- matrix(NA, nrow=nsnp(data), ncol=nbatch, dimnames=list(snpID, exp.batch))
  for (i in 1:nsnp(data)) {
    for (j in 1:nbatch) {
      tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                      nrow=2, ncol=2)
      try({
        chisq <- chisq.test(tbl, correct=TRUE)
        exp.chisq[i,j] <- chisq[["statistic"]]
      })
    }
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- rep(NA, nbatch)
  names(exp.ave) <- exp.batch
  exp.lam <- rep(NA, nbatch)
  names(exp.lam) <- exp.batch
  for (i in 1:nbatch) {
    exp.ave[i] <- mean(exp.chisq[,i], na.rm=TRUE)
    exp.lam[i] <- median(exp.chisq[,i], na.rm=TRUE) / 0.456
  }
  
  # test function
  resfile <- tempfile()
  batchChisqTest(data, batchVar="batch", return.by.snp=TRUE, outfile=resfile)
  res <- getobj(paste(resfile, "RData", sep="."))
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)
  unlink(paste(resfile, "*", sep=""))
}

checkExclude <- function(nc) {
  # define batches - 4
  scanID <- 1:20
  batch <- rep(paste("batch", 1:4, sep=""), 5)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)
  scan.exclude <- c(1,2,10)

  # expected values
  snpID <- getSnpID(data)
  geno <- getGenotype(data)[, !(scanID %in% scan.exclude)]
  batch <- batch[!(scanID %in% scan.exclude)]
  exp.batch <- sort(unique(batch))
  nbatch <- length(exp.batch)
  nA <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  nB <- matrix(0, nrow=nsnp(data), ncol=nbatch)
  for (i in 1:nbatch) {
    thisbatch <- batch == exp.batch[i]
    nA[,i] <- rowSums(geno[,thisbatch], na.rm=TRUE)
    nB[,i] <- rowSums(2-geno[,thisbatch], na.rm=TRUE)
  }
  exp.chisq <- matrix(NA, nrow=nsnp(data), ncol=nbatch, dimnames=list(snpID, exp.batch))
  for (i in 1:nsnp(data)) {
    for (j in 1:nbatch) {
      tbl <- matrix(c(nA[i,j], nB[i,j], sum(nA[i,-j]), sum(nB[i,-j])),
                      nrow=2, ncol=2)
      try({
        chisq <- chisq.test(tbl, correct=TRUE)
        exp.chisq[i,j] <- chisq[["statistic"]]
      })
    }
  }
  exp.chisq[!is.finite(exp.chisq)] <- NA
  exp.chisq[pmin(rowSums(nA), rowSums(nB)) == 0] <- NA
  exp.ave <- rep(NA, nbatch)
  names(exp.ave) <- exp.batch
  exp.lam <- rep(NA, nbatch)
  names(exp.lam) <- exp.batch
  for (i in 1:nbatch) {
    exp.ave[i] <- mean(exp.chisq[,i], na.rm=TRUE)
    exp.lam[i] <- median(exp.chisq[,i], na.rm=TRUE) / 0.456
  }

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=TRUE,
                        scan.exclude=scan.exclude)
  checkEquals(exp.ave, res$mean.chisq)
  checkEquals(exp.lam, res$lambda)
  checkEquals(exp.chisq, res$chisq, tolerance=1e-7)

}

checkNoSnp <- function(nc) {
  # define batches - 4
  scanID <- 1:20
  batch <- rep(paste("batch", 1:4, sep=""), 5)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex, batch=batch)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  data <- GenotypeData(nc, scanAnnot=scanAnnot)

  # test function
  res <- batchChisqTest(data, batchVar="batch", return.by.snp=FALSE)
  checkEquals(2, length(res))
  checkIdentical(c("mean.chisq", "lambda"), names(res))
}


test_batchChisqTest <- function() {
  # simulate data - autosomes only
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=22,
                           n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)

  # run tests
  check2batches(nc)
  check4batches(nc)
  checkNoYates(nc)
  checkFileOut(nc)
  checkExclude(nc)
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
