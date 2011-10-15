test_ncdfAddIntensity <- function() {
  # first create empty netCDF
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  ncdfCreate(snpAnnot, ncfile, variables=c("quality","X","Y"),
                  n.samples=nrow(scanAnnot))

  # add sampleID and quality
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")
  scanAnnot1 <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot1) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(2,4)); names(col.nums) <- c("snp", "qs")
  diagfile <- tempfile()
  res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot1, sep.type="\t",
                       skip.num=1, col.total=6, col.nums=col.nums,
                       scan.name.in.file=-1, diagnostics.filename=diagfile)
  file.remove(diagfile)

  # add intensity
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  res <- ncdfAddIntensity(path, ncfile, snpAnnot, scanAnnot)
  checkTrue(all(res$chk == 1))

  # check
  nc <- NcdfIntensityReader(ncfile)
  origfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  nc2 <- NcdfIntensityReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkEquals(getX(nc), getX(nc2, snp=c(1,-1), scan=c(1,3)), tolerance=1e-7)
  checkEquals(getY(nc), getY(nc2, snp=c(1,-1), scan=c(1,3)), tolerance=1e-7)
  close(nc)
  close(nc2)
  
  # check error conditions 
  names(snpAnnot) <- c("snpID", "snpName")
  names(scanAnnot) <- c("scanID", "foo", "file") # bad scan names
  checkException({
    res <- ncdfAddIntensity(path, ncfile, snpAnnot, scanAnnot)
  })
  
  names(scanAnnot) <- c("scanID", "scanName", "file")
  snpAnnot$snpID[1] <- snpAnnot$snpID[2] # bad snp IDs
  checkException({
    res <- ncdfAddIntensity(path, ncfile, snpAnnot, scanAnnot)
  })
  
  file.remove(ncfile)
}
