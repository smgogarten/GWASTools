test_qualityScoreByScan <- function() {
  # simulate data
  intenfile <- tempfile()
  simulateIntensityMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=intenfile)
  intennc <- GdsIntensityReader(intenfile)
  scanID <- 1:20
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  intenData <- IntensityData(intennc, scanAnnot=scanAnnot)
  
  genofile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=genofile)
  genonc <- GdsGenotypeReader(genofile)
  genoData <- GenotypeData(genonc, scanAnnot=scanAnnot)

  # expected results
  qual <- getQuality(intenData)
  geno <- getGenotype(genoData)
  snpID <- getSnpID(genoData)
  chrom <- getChromosome(genoData, char=TRUE)
  qual[is.na(geno)] <- NA
  qual[chrom == "Y", sex == "F"] <- NA
  mean.quality <- colMeans(qual, na.rm=TRUE)
  median.quality <- rowMedians(t(qual), na.rm=TRUE)
  exp <- cbind(mean.quality, median.quality)
  rownames(exp) <- scanID

  # test
  res <- qualityScoreByScan(intenData, genoData) 
  checkEquals(exp, res)
  
  # snp.exclude - expected results
  snp.exclude <- c(1,2,10)
  qual <- qual[c(-1,-2,-10),]
  
  mean.quality <- colMeans(qual, na.rm=TRUE)
  median.quality <- rowMedians(t(qual), na.rm=TRUE)
  exp <- cbind(mean.quality, median.quality)
  rownames(exp) <- scanID
  
  # test
  res <- qualityScoreByScan(intenData, genoData, snp.exclude=snp.exclude) 
  checkEquals(exp, res)

  close(intenData)
  close(genoData)
  file.remove(intenfile, genofile)
}
