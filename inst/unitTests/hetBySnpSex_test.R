test_hetBySnpSex <- function() {  
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

  # expected results
  geno <- getGenotype(genoData)
  snpID <- getSnpID(genoData)
  males <- sex == "M"
  het.M <- rowSums((geno[,males] == 1), na.rm=TRUE) / rowSums(!is.na(geno[,males]))
  names(het.M) <- snpID
  females <- sex == "F"
  het.F <- rowSums((geno[,females] == 1), na.rm=TRUE) / rowSums(!is.na(geno[,females]))
  names(het.F) <- snpID

  # test
  het <- hetBySnpSex(genoData)
  checkEquals(het[,"M"], het.M)
  checkEquals(het[,"F"], het.F)
  
  # scan.exclude - expected results
  scan.exclude <- c(1,2,10)
  males <- sex == "M" & !(scanID %in% scan.exclude)
  het.M <- rowSums((geno[,males] == 1), na.rm=TRUE) / rowSums(!is.na(geno[,males]))
  names(het.M) <- snpID
  females <- sex == "F" & !(scanID %in% scan.exclude)
  het.F <- rowSums((geno[,females] == 1), na.rm=TRUE) / rowSums(!is.na(geno[,females]))
  names(het.F) <- snpID

  # test
  het <- hetBySnpSex(genoData, scan.exclude=scan.exclude)
  checkEquals(het[,"M"], het.M)
  checkEquals(het[,"F"], het.F)
  
  close(genoData)
  file.remove(ncfile)
}
