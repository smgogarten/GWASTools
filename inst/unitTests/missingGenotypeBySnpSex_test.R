test_missingGenotypeBySnpSex <- function() {
  # simulate data
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=ncfile)
  nc <- GdsGenotypeReader(ncfile)
  scanID <- 1:20
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot)
  ychr <- getChromosome(genoData, char=TRUE) == "Y"

  # expected results
  geno <- getGenotype(genoData)
  snpID <- getSnpID(genoData)
  males <- sex == "M"
  miss.M <- rowSums(is.na(geno[,males]))
  names(miss.M) <- snpID
  females <- sex == "F"
  miss.F <- rowSums(is.na(geno[,females]))
  names(miss.F) <- snpID
  miss.frac <- rowSums(is.na(geno)) / length(scanID)
  miss.frac[ychr] <- miss.M[ychr] / sum(males)
  names(miss.frac) <- snpID

  # test
  miss <- missingGenotypeBySnpSex(genoData)
  checkEquals(miss$missing.counts[,"M"], miss.M)
  checkEquals(miss$missing.counts[,"F"], miss.F)
  checkEquals(miss$scans.per.sex["M"], sum(males), checkNames=FALSE)
  checkEquals(miss$scans.per.sex["F"], sum(females), checkNames=FALSE)
  checkEquals(miss$missing.fraction, miss.frac)

  # scan.exclude - expected results
  scan.exclude <- c(1,2,10)
  males <- sex == "M" & !(scanID %in% scan.exclude)
  miss.M <- rowSums(is.na(geno[,males]))
  names(miss.M) <- snpID
  females <- sex == "F" & !(scanID %in% scan.exclude)
  miss.F <- rowSums(is.na(geno[,females]))
  names(miss.F) <- snpID
  miss.frac <- rowSums(is.na(geno[,(males | females)])) /
    (length(scanID) - length(scan.exclude))
  miss.frac[ychr] <- miss.M[ychr] / sum(males)
  names(miss.frac) <- snpID

  # test
  miss <- missingGenotypeBySnpSex(genoData, scan.exclude=scan.exclude)
  checkEquals(miss$missing.counts[,"M"], miss.M)
  checkEquals(miss$missing.counts[,"F"], miss.F)
  checkEquals(miss$scans.per.sex["M"], sum(males), checkNames=FALSE)
  checkEquals(miss$scans.per.sex["F"], sum(females), checkNames=FALSE)
  checkEquals(miss$missing.fraction, miss.frac)

  close(genoData)
  file.remove(ncfile)
}
