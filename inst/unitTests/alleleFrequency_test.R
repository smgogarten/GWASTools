test_alleleFrequency <- function() {
  # simulate data
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)
  scanID <- 1:20
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot)
  chrom <- getChromosome(genoData, char=TRUE)
  males <- sex == "M"
  females <- sex == "F"

  # expected results - autosomes only
  auto <- chrom %in% c(1:22)
  geno <- getGenotype(genoData)[auto,]
  snpID <- getSnpID(genoData)[auto]
  afreq.M <- rowSums(geno[,males], na.rm=TRUE) /
    (2*rowSums(!is.na(geno[,males])))
  names(afreq.M) <- snpID
  afreq.F <- rowSums(geno[,females], na.rm=TRUE) /
    (2*rowSums(!is.na(geno[,females])))
  names(afreq.F) <- snpID
  afreq.all <- rowSums(geno, na.rm=TRUE) /
    (2*rowSums(!is.na(geno)))
  names(afreq.all) <- snpID

  # test
  afreq <- alleleFrequency(genoData)
  checkEquals(afreq[auto,"M"], afreq.M)
  checkEquals(afreq[auto,"F"], afreq.F)
  checkEquals(afreq[auto,"all"], afreq.all)
  checkEquals(afreq[auto,"MAF"], pmin(afreq.all, 1-afreq.all))
  checkTrue(max(afreq[auto,"MAF"], na.rm=TRUE) <= 0.5)
  
  # X chrom - males
  x <- chrom == "X"
  geno <- getGenotype(genoData)[x, males]
  snpID <- getSnpID(genoData)[x]
  geno[geno == 1] <- NA # don't count hets
  acnt.M <- rowSums(geno, na.rm=TRUE) / 2
  nonmiss.M <-  rowSums(!is.na(geno))
  afreq.M <- acnt.M / nonmiss.M
  names(afreq.M) <- snpID

  # X chrom - females
  geno <- getGenotype(genoData)[x, females]
  acnt.F <- rowSums(geno, na.rm=TRUE)
  nonmiss.F <- 2*rowSums(!is.na(geno))
  afreq.F <-  acnt.F / nonmiss.F
  names(afreq.F) <- snpID

  # X chrom - all
  afreq.all <- (acnt.M + acnt.F) / (nonmiss.M + nonmiss.F)
  names(afreq.all) <- snpID
  
  # test
  checkEquals(afreq[x,"M"], afreq.M)
  checkEquals(afreq[x,"F"], afreq.F)
  checkEquals(afreq[x,"all"], afreq.all)
  checkEquals(afreq[x,"MAF"], pmin(afreq.all, 1-afreq.all))
  checkTrue(max(afreq[x,"MAF"], na.rm=TRUE) <= 0.5)
  
  # Y chrom - males
  y <- chrom == "Y"
  geno <- getGenotype(genoData)[y, males]
  snpID <- getSnpID(genoData)[y]
  geno[geno == 1] <- NA # don't count hets
  acnt.M <- rowSums(geno, na.rm=TRUE) / 2
  nonmiss.M <-  rowSums(!is.na(geno))
  afreq.M <- acnt.M / nonmiss.M
  names(afreq.M) <- snpID

  # Y chrom - females
  geno <- getGenotype(genoData)[y, females]
  geno[,] <- NA
  acnt.F <- rowSums(geno, na.rm=TRUE)
  nonmiss.F <- 2*rowSums(!is.na(geno))
  afreq.F <-  acnt.F / nonmiss.F
  names(afreq.F) <- snpID

  # Y chrom - all
  afreq.all <- (acnt.M + acnt.F) / (nonmiss.M + nonmiss.F)
  names(afreq.all) <- snpID
  
  # test
  checkEquals(afreq[y,"M"], afreq.M)
  checkEquals(afreq[y,"F"], afreq.F)
  checkEquals(afreq[y,"all"], afreq.all)
  checkEquals(afreq[y,"MAF"], pmin(afreq.all, 1-afreq.all))
  checkTrue(max(afreq[y,"MAF"], na.rm=TRUE) <= 0.5)
  
  # scan.exclude - expected results
  scan.exclude <- c(1,2,10)
  males <- sex == "M" & !(scanID %in% scan.exclude)
  females <- sex == "F" & !(scanID %in% scan.exclude)
  all <- !(scanID %in% scan.exclude)
  auto <- chrom %in% c(1:22)
  geno <- getGenotype(genoData)[auto,]
  snpID <- getSnpID(genoData)[auto]
  afreq.M <- rowSums(geno[,males], na.rm=TRUE) /
    (2*rowSums(!is.na(geno[,males])))
  names(afreq.M) <- snpID
  afreq.F <- rowSums(geno[,females], na.rm=TRUE) /
    (2*rowSums(!is.na(geno[,females])))
  names(afreq.F) <- snpID
  afreq.all <- rowSums(geno[,all], na.rm=TRUE) /
    (2*rowSums(!is.na(geno[,all])))
  names(afreq.all) <- snpID

  # test
  afreq <- alleleFrequency(genoData, scan.exclude=scan.exclude)
  checkEquals(afreq[auto,"M"], afreq.M)
  checkEquals(afreq[auto,"F"], afreq.F)
  checkEquals(afreq[auto,"all"], afreq.all)
  checkEquals(afreq[auto,"MAF"], pmin(afreq.all, 1-afreq.all))
  checkTrue(max(afreq[auto,"MAF"], na.rm=TRUE) <= 0.5)
  
  close(genoData)
  file.remove(ncfile)
}


.testData <- function(chromosome, nsamp=100) {
    nsnp <- length(chromosome)
    geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE), nrow=nsnp, ncol=nsamp)
    geno[1,] <- 1 # at least one snp all non-missing
    mgr <- MatrixGenotypeReader(geno, snpID=1:nsnp, scanID=1:nsamp,
                                chromosome=as.integer(chromosome),
                                position=1:nsnp)

    scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=1:nsamp,
                                         sex=sample(c("M","F",NA), nsamp, replace=TRUE),
                                         trait=rnorm(nsamp, mean=10, sd=2),
                                         stringsAsFactors=FALSE))

    GenotypeData(mgr, scanAnnot=scanAnnot)
}

test_missingSex <- function() {
    chr <- c(rep(1,50), rep(23,30), rep(25,20))
    genoData <- .testData(chr)
    sex <- getSex(genoData)
    male <- sex %in% "M"
    female <- sex %in% "F"
    afreq <- alleleFrequency(genoData, verbose=FALSE)
    checkTrue(all(afreq[,"n.M"] + afreq[,"n.F"] <= afreq[,"n"]))
    checkEquals(sum(male), max(afreq[,"n.M"]))
    checkEquals(sum(female), max(afreq[,"n.F"]))
    checkEquals(nscan(genoData), max(afreq[,"n"]))

    geno <- getGenotype(genoData)
    auto <- chr %in% 1:22
    xchr <- chr %in% 23
    ychr <- chr %in% 25
    for (i in which(xchr | ychr)) geno[i,male][geno[i,male] %in% 1] <- NA
    checkEquals(afreq[auto,"M"], 0.5*rowMeans(geno[auto,male], na.rm=TRUE), checkNames=FALSE)
    checkEquals(afreq[auto,"F"], 0.5*rowMeans(geno[auto,female], na.rm=TRUE), checkNames=FALSE)
    checkEquals(afreq[auto,"all"], 0.5*rowMeans(geno[auto,], na.rm=TRUE), checkNames=FALSE)
    checkEquals(afreq[xchr,"F"], 0.5*rowMeans(geno[xchr,female], na.rm=TRUE), checkNames=FALSE)
    checkEquals(afreq[xchr,"M"], 0.5*rowMeans(geno[xchr,male], na.rm=TRUE), checkNames=FALSE)
    checkEquals(afreq[xchr,"all"], (0.5*rowSums(geno[xchr,male], na.rm=TRUE) + rowSums(geno[xchr,female], na.rm=TRUE)) /
                (afreq[xchr,"n.M"] + 2*afreq[xchr,"n.F"]), checkNames=FALSE)
    checkTrue(all(is.na(afreq[ychr,"F"])))
    checkEquals(afreq[ychr,"M"], 0.5*rowMeans(geno[ychr,male], na.rm=TRUE), checkNames=FALSE)
    checkEquals(afreq[ychr,"all"], afreq[ychr,"M"], checkNames=FALSE)
    
}

test_compare <- function() {
    # autosomal - everything should match
    genoData <- .testData(rep(1,100))
    afreq <- alleleFrequency(genoData, verbose=FALSE)
    hwe <- exactHWE(genoData, verbose=FALSE)
    assoc <- assocRegression(genoData, "trait", verbose=FALSE)
    checkEquals(afreq[,"MAF"], hwe[,"MAF"], checkNames=FALSE)
    checkEquals(afreq[,"MAF"], assoc[,"MAF"], checkNames=FALSE)

    # xchrom - females only for HWE
    genoData <- .testData(rep(23,100))
    afreq <- alleleFrequency(genoData, verbose=FALSE)
    hwe <- exactHWE(genoData, verbose=FALSE)
    checkEquals(pmin(afreq[,"F"], 1-afreq[,"F"]), hwe[,"MAF"], checkNames=FALSE)
}
