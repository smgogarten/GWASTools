.makeTestFile <- function() {
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  data(illumina_snp_annot, illumina_scan_annot)
  snp <- illumina_snp_annot
  samp <- illumina_scan_annot
  samp <- samp[samp$file %in% list.files(path),]
  gdsfile <- tempfile()

  snp <- snp[,c("snpID", "rsID", "chromosome", "position", "alleleA", "alleleB")]
  names(snp)[1:2] <- c("snpID", "snpName")
  samp <- samp[,c("scanID", "genoRunID", "file")]
  names(samp) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,5,16,17))
  names(col.nums) <- c("snp", "sample", "quality", "X", "Y")
  diagfile <- tempfile()
  res <- createDataFile(path, gdsfile, file.type="gds", variables=c("quality","X","Y"),
                        snp, samp, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  unlink(diagfile)
  return(gdsfile)
}


test_GdsIntensityReader <- function() {
  gdsfile <- .makeTestFile()
  obj <- GdsIntensityReader(gdsfile)

  x <- getX(obj)
  checkIdentical(c(nsnp(obj),nscan(obj)), dim(x))

  nsnp <- 100L
  nsamp <- 2L
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
  sel <- scanID %in% sample(scanID, 2)
  checkIdentical(scanID[sel], getScanID(obj, sel))

  chromChar <- getChromosome(obj, char=TRUE)
  checkTrue(is.character(chromChar))
  checkTrue(all(chromChar %in% c(1:22,"X","Y","XY","M","U")))
  close(obj)

  unlink(gdsfile)

  ## # check exception with incorrect dimensions
  ## dim1 <- dim.def.ncdf("dim1", "count", 1:10)
  ## var1 <- var.def.ncdf("var1", "count", dim=dim1, missval=-1)
  ## file <- tempfile()
  ## nc <- create.ncdf(file, var1)
  ## close.ncdf(nc)
  ## checkException(GdsIntensityReader(file))
  ## unlink(file)

  ## # check exception with incorrect variable names
  ## snp <- dim.def.ncdf("snp", "count", 1:10)
  ## sampleID <- var.def.ncdf("sampleID", "count", dim=snp, missval=-1)
  ## file <- tempfile()
  ## nc <- create.ncdf(file, sampleID)
  ## close.ncdf(nc)
  ## checkException(GdsIntensityReader(file))
  ## unlink(file)
}

test_GdsIntensityReader_error <- function() {
  
  gdsfile <- .makeTestFile()
  # make sure it can be opened
  obj <- GdsIntensityReader(gdsfile)
  close(obj)
  
  checkException(GdsIntensityReader(gdsfile, positionVar="pos"))
  # make sure it can be opened
  obj <- GdsIntensityReader(gdsfile)
  close(obj)
  
  # now check that it's not closed when there is an error and a gds object is passed
  gds <- openfn.gds(gdsfile)
  checkException(GdsIntensityReader(gds, positionVar="pos"))
  obj <- GdsIntensityReader(gds)
  
  close(obj)
  unlink(gdsfile)
}
