test_NcdfIntensityReader <- function() {  
  file <- tempfile()
  simulateIntensityMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=file, file.type="ncdf")
  obj <- NcdfIntensityReader(file)
  
  x <- getX(obj)
  checkIdentical(c(nsnp(obj),nscan(obj)), dim(x))
  
  nsnp <- 100L
  nsamp <- 10L
  x <- getX(obj, snp=c(1,nsnp), scan=c(1,nsamp))
  checkIdentical(c(nsnp,nsamp), dim(x))
  y <- getY(obj, snp=c(1,nsnp), scan=c(1,nsamp))
  checkIdentical(c(nsnp,nsamp), dim(y))
  q <- getQuality(obj, snp=c(1,nsnp), scan=c(1,nsamp))
  checkIdentical(c(nsnp,nsamp), dim(q))

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
  obj <- NcdfIntensityReader(file, YchromCode=24L, XYchromCode=25L)
  chromChar <- getChromosome(obj, char=TRUE)
  checkIdentical(rep(c(1:22,"X","Y","XY","M"), each=10), chromChar)

  close(obj)
  unlink(file)

  file <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
  obj <- NcdfIntensityReader(file)
  
  nsnp <- 100L
  nsamp <- 10L
  baf <- getBAlleleFreq(obj, snp=c(1,nsnp), scan=c(1,nsamp))
  checkIdentical(c(nsnp,nsamp), dim(baf))
  lrr <- getLogRRatio(obj, snp=c(1,nsnp), scan=c(1,nsamp))
  checkIdentical(c(nsnp,nsamp), dim(lrr))

  close(obj)  

  # check exception with incorrect dimensions
  dim1 <- ncdim_def("dim1", "count", 1:10)
  var1 <- ncvar_def("var1", "count", dim=dim1, missval=-1)
  file <- tempfile()
  nc <- nc_create(file, var1)
  nc_close(nc)
  checkException(NcdfIntensityReader(file))
  unlink(file)
  
  # check exception with incorrect variable names
  snp <- ncdim_def("snp", "count", 1:10)
  sampleID <- ncvar_def("sampleID", "count", dim=snp, missval=-1)
  file <- tempfile()
  nc <- nc_create(file, sampleID)
  nc_close(nc)
  checkException(NcdfIntensityReader(file))
  unlink(file)
}
