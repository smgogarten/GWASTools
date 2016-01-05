test_hetByScanChrom <- function() {
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
  nchr <- 26
  het.cnt <- matrix(NA, nrow=length(scanID), ncol=(nchr+1),
                    dimnames=list(scanID, c(uchr,"A")))
  nm.cnt <- matrix(NA, nrow=length(scanID), ncol=(nchr+1),
                     dimnames=list(scanID, c(uchr,"A")))
  for (i in 1:nchr) {
    gc <- geno[(chrom == uchr[i]),]
    het.cnt[,i] <- colSums(gc == 1, na.rm=TRUE)
    nm.cnt[,i] <- colSums(!is.na(gc))
  }
  ga <- geno[(chrom %in% 1:22),]
  het.cnt[,"A"] <- colSums(ga == 1, na.rm=TRUE)
  nm.cnt[,"A"] <- colSums(!is.na(ga))
  het.frac <- het.cnt / nm.cnt

  het <- hetByScanChrom(genoData)
  checkEquals(het, het.frac)

  # snp.exclude - expected results
  snpID <- getSnpID(genoData)
  snp.exclude <- snpID[c(1,2,11)]
  het.cnt <- matrix(NA, nrow=length(scanID), ncol=(nchr+1),
                    dimnames=list(scanID, c(uchr,"A")))
  nm.cnt <- matrix(NA, nrow=length(scanID), ncol=(nchr+1),
                     dimnames=list(scanID, c(uchr,"A")))
  for (i in 1:nchr) {
    gc <- geno[(chrom == uchr[i] & !(snpID %in% snp.exclude)),]
    het.cnt[,i] <- colSums(gc == 1, na.rm=TRUE)
    nm.cnt[,i] <- colSums(!is.na(gc))
  }
  ga <- geno[(chrom %in% 1:22 & !(snpID %in% snp.exclude)),]
  het.cnt[,"A"] <- colSums(ga == 1, na.rm=TRUE)
  nm.cnt[,"A"] <- colSums(!is.na(ga))
  het.frac <- het.cnt / nm.cnt
  
  het <- hetByScanChrom(genoData, snp.exclude=snp.exclude)
  checkEquals(het, het.frac)
  
  close(genoData)
  file.remove(ncfile)
}
