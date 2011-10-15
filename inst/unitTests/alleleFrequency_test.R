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
  
  close(genoData)
  file.remove(ncfile)
}
