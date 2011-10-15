test_qualityScoreBySnp <- function() {
  # simulate data
  intenfile <- tempfile()
  simulateIntensityMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=intenfile)
  intennc <- NcdfIntensityReader(intenfile)
  scanID <- 1:20
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- data.frame(scanID=scanID, sex=sex)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  intenData <- IntensityData(intennc, scanAnnot=scanAnnot)
  
  genofile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=genofile)
  genonc <- NcdfGenotypeReader(genofile)
  genoData <- GenotypeData(genonc, scanAnnot=scanAnnot)

  # expected results
  qual <- getQuality(intenData)
  geno <- getGenotype(genoData)
  snpID <- getSnpID(genoData)
  qual[is.na(geno)] <- NA
  mean.quality <- rowMeans(qual, na.rm=TRUE)
  median.quality <- rowMedians(qual, na.rm=TRUE)
  exp <- cbind(mean.quality, median.quality)
  rownames(exp) <- snpID

  # test
  res <- qualityScoreBySnp(intenData, genoData) 
  checkEquals(exp, res)
  
  # scan.exclude - expected results
  scan.exclude <- c(1,2,10)
  qual <- qual[,c(-1,-2,-10)]
  
  mean.quality <- rowMeans(qual, na.rm=TRUE)
  median.quality <- rowMedians(qual, na.rm=TRUE)
  exp <- cbind(mean.quality, median.quality)
  rownames(exp) <- snpID

  # test
  res <- qualityScoreBySnp(intenData, genoData, scan.exclude=scan.exclude) 
  checkEquals(exp, res)

  # test block size
  res <- qualityScoreBySnp(intenData, genoData, scan.exclude=scan.exclude,
                           block.size=50) 
  checkEquals(exp, res)
  
  # test block size - 1 snp at a time
  res <- qualityScoreBySnp(intenData, genoData, scan.exclude=scan.exclude,
                           block.size=1) 
  checkEquals(exp, res)
  
  close(intenData)
  close(genoData)
  file.remove(intenfile, genofile)
}
