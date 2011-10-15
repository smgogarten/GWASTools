test_NcdfGenotypeReader <- function() {  
  file <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=file)
  obj <- NcdfGenotypeReader(file)
  
  geno <- getGenotype(obj)
  checkIdentical(c(nsnp(obj),nscan(obj)), dim(geno))
  
  nsnp <- 100L
  nsamp <- 10L
  geno <- getGenotype(obj, snp=c(1,nsnp), scan=c(1,nsamp))
  checkIdentical(c(nsnp,nsamp), dim(geno))

  checkIdentical(length(getSnpID(obj)), nsnp(obj))
  checkIdentical(length(getChromosome(obj)), nsnp(obj))
  checkIdentical(length(getPosition(obj)), nsnp(obj))
  checkIdentical(length(getScanID(obj)), nscan(obj))

  scanID <- getScanID(obj)
  sel <- scanID %in% sample(scanID, 10)
  checkIdentical(scanID[sel], getScanID(obj, sel))

  chromChar <- getChromosome(obj, char=TRUE)
  checkTrue(is.character(chromChar))
  checkTrue(all(chromChar %in% c(1:22,"X","Y","XY","M","U")))
  checkIdentical(rep(c(1:22,"X","XY","Y","M"), each=10), chromChar)
  close(obj)

  # check alternate chromosome codes
  obj <- NcdfGenotypeReader(file, YchromCode=24L, XYchromCode=25L)
  chromChar <- getChromosome(obj, char=TRUE)
  checkIdentical(rep(c(1:22,"X","Y","XY","M"), each=10), chromChar)
  close(obj)
  unlink(file)

  # check exception with incorrect dimensions
  dim1 <- dim.def.ncdf("dim1", "count", 1:10)
  var1 <- var.def.ncdf("var1", "count", dim=dim1, missval=-1)
  file <- tempfile()
  nc <- create.ncdf(file, var1)
  close.ncdf(nc)
  checkException(NcdfGenotypeReader(file))
  unlink(file)
  
  # check exception with incorrect variable names
  snp <- dim.def.ncdf("snp", "count", 1:10)
  sampleID <- var.def.ncdf("sampleID", "count", dim=snp, missval=-1)
  file <- tempfile()
  nc <- create.ncdf(file, sampleID)
  close.ncdf(nc)
  checkException(NcdfGenotypeReader(file))
  unlink(file)
}
