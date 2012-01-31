test_plinkCheck <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc.orig <- NcdfGenotypeReader(ncfile)
  
  scanID <- getScanID(nc.orig)
  sex <- c(rep("M", 10), rep("F", 10))
  family <- rep(0, 20)
  father <- rep(0, 20)
  mother <- rep(0, 20)
  scandf <- data.frame(family, scanID, father, mother, sex, stringsAsFactors=FALSE)
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  rsID <- paste("rs", snpID, sep="")
  allele.A <- sample(c("A","C","G","T"), nsnp(nc.orig), replace=TRUE)
  allele.B <- sample(c("A","C","G","T"), nsnp(nc.orig), replace=TRUE)
  snpdf <- data.frame(snpID, chromosome, position, rsID, allele.A,
                      allele.B, stringsAsFactors=FALSE)
  geno <- getGenotype(nc.orig)
  genoData <- GenotypeData(nc.orig, scanAnnot=ScanAnnotationDataFrame(scandf),
                           snpAnnot=SnpAnnotationDataFrame(snpdf))
  
  pedfile <- tempfile()
  plinkWrite(genoData, pedfile, alleleA.col="allele.A", alleleB.col="allele.B")
  logfile <- tempfile()
  checkTrue(plinkCheck(genoData, pedfile, logfile, alleleA.col="allele.A", alleleB.col="allele.B"))
  close(genoData)

  # change ncdf file
  nc.new <- open.ncdf(ncfile, write=TRUE)
  geno <- get.var.ncdf(nc.new, "genotype")
  geno[geno == 1] <- 2
  put.var.ncdf(nc.new, "genotype", geno)
  close.ncdf(nc.new)
  nc.new <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc.new, scanAnnot=ScanAnnotationDataFrame(scandf),
                           snpAnnot=SnpAnnotationDataFrame(snpdf))
  checkTrue(!plinkCheck(genoData, pedfile, logfile, alleleA.col="allele.A", alleleB.col="allele.B"))
  close(genoData)
  
  unlink(c(ncfile, logfile, paste(pedfile, "*", sep=".")))
}
