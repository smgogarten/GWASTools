test_ncdf <- function() {
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "probeID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  diagfile <- tempfile()
  res <- createAffyIntensityFile(path, ncfile, file.type="ncdf",
                                 snpAnnot, scanAnnot,
                                 diagnostics.filename=diagfile)

  # check
  nc <- NcdfIntensityReader(ncfile)
  origfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  nc2 <- NcdfIntensityReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkIdentical(getX(nc), getX(nc2, snp=c(1,-1), scan=c(1,3)))
  checkIdentical(getY(nc), getY(nc2, snp=c(1,-1), scan=c(1,3)))
  close(nc)
  close(nc2)

  file.remove(diagfile)
  file.remove(ncfile)
}


test_gds <- function() {
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "probeID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  diagfile <- tempfile()
  res <- createAffyIntensityFile(path, ncfile, file.type="gds",
                                 snpAnnot, scanAnnot,
                                 diagnostics.filename=diagfile)

  # check
  nc <- GdsIntensityReader(ncfile)
  origfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  nc2 <- NcdfIntensityReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkIdentical(getX(nc), getX(nc2, snp=c(1,-1), scan=c(1,3)))
  checkIdentical(getY(nc), getY(nc2, snp=c(1,-1), scan=c(1,3)))
  close(nc)
  close(nc2)

  file.remove(diagfile)
  file.remove(ncfile)
}
