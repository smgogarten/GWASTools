test_asSnpMatrix <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=2,
                         n.samples=5, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)
  scanID <- getScanID(nc)
  subjID <- paste("A", scanID, sep="")
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  snpID <- getSnpID(nc)
  chrom <- getChromosome(nc)
  pos <- getPosition(nc)
  rsID <- paste("rs", snpID, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      rsID=rsID, stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  geno <- getGenotype(genoData)
  gc <- matrix("NA", nrow=nsnp(genoData), ncol=nscan(genoData))
  gc[geno %in% 0] <- "B/B"
  gc[geno %in% 1] <- "A/B"
  gc[geno %in% 2] <- "A/A"

  snpmat <- asSnpMatrix(genoData)
  gcn <- t(gc)
  dimnames(gcn) <- list(scanID, snpID)
  checkIdentical(gcn, as(snpmat, "character"))

  # use different names for scan and snp
  snpmat <- asSnpMatrix(genoData, snpNames="rsID", scanNames="subjID")
  gcn <- t(gc)
  dimnames(gcn) <- list(subjID, rsID)
  checkIdentical(gcn, as(snpmat, "character"))

  # use a subset of data
  snpmat <- asSnpMatrix(genoData, snp=c(2,10), scan=c(3,-1))
  gcn <- t(gc[2:11, 3:5])
  dimnames(gcn) <- list(scanID[3:5], snpID[2:11])
  checkIdentical(gcn, as(snpmat, "character"))

  close(genoData)
  unlink(ncfile)
}
