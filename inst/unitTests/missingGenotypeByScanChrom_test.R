test_missingGenotypeByScanChrom <- function() {
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
  chrom <- getChromosome(genoData, char=TRUE)
  uchr <- unique(chrom)
  geno <- getGenotype(genoData)

  # expected results
  spc <- rep(10, length(uchr))
  names(spc) <- uchr
  miss.cnt <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  for (i in 1:length(uchr)) {
    miss.cnt[,i] <- colSums(is.na(geno[(chrom == uchr[i]),]))
  }
  miss.frac <- colSums(is.na(geno)) / length(chrom)
  female <- sex == "F"
  notY <- chrom != "Y"
  miss.frac[female] <- colSums(is.na(geno)[notY, female]) /
    sum(notY)
  names(miss.frac) <- scanID

  miss <- missingGenotypeByScanChrom(genoData)
  checkEquals(miss$missing.counts, miss.cnt)
  checkEquals(miss$snps.per.chr, spc)
  checkEquals(miss$missing.fraction, miss.frac)

  # snp.exclude - expected results
  snpID <- getSnpID(genoData)
  snp.exclude <- snpID[c(1,2,11)]
  spc <- c(8,9,rep(10,24))
  names(spc) <- uchr
  miss.cnt <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  for (i in 1:length(uchr)) {
    miss.cnt[,i] <- colSums(is.na(geno[(chrom == uchr[i] &
                                        !(snpID %in% snp.exclude)),]))
  }
  miss.frac <- colSums(is.na(geno[!(snpID %in% snp.exclude),])) /
    (length(chrom) - length(snp.exclude))
  female <- sex == "F"
  notY <- chrom != "Y" & !(snpID %in% snp.exclude)
  miss.frac[female] <- colSums(is.na(geno)[notY, female]) / sum(notY)
  names(miss.frac) <- scanID
  
  miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snp.exclude)
  checkEquals(miss$missing.counts, miss.cnt)
  checkEquals(miss$snps.per.chr, spc)
  checkEquals(miss$missing.fraction, miss.frac)
  
  close(genoData)
  file.remove(ncfile)
}
