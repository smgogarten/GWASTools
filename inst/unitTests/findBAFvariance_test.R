test_findBAFvariance <- function() {
  blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
  blnc <- NcdfIntensityReader(blfile)

  genofile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
  genonc <- NcdfGenotypeReader(genofile)

  nbins <- rep(8, 3) # chroms 21-26 in this file, so need bins for (21,22,23)
  baf.res <- sdByScanChromWindow(blnc, genonc, nbins=nbins)
  checkEquals(length(baf.res), 3)
  checkEquals(dim(baf.res[[1]]), c(nscan(blnc), 7))

  data(illumina_scan_annot)
  sex <- illumina_scan_annot$sex
  sd.res <- meanSdByChromWindow(baf.res, sex)

  var.res <- findBAFvariance(sd.res, baf.res, sex, 2)

  sd.res <- medianSdOverAutosomes(baf.res)
  checkEquals(dim(sd.res), c(nscan(blnc), 2))

  # default value for nbins
  baf.res <-  sdByScanChromWindow(blnc, genonc)
  checkEquals(dim(baf.res[[1]]), c(nscan(blnc), 1))

  sd.res <- medianSdOverAutosomes(baf.res)
  checkEquals(dim(sd.res), c(nscan(blnc), 2))

  # try LRR
  lrr.res <-  sdByScanChromWindow(blnc, var="LogRRatio", incl.hom=TRUE)
  checkEquals(length(lrr.res), 3)
  checkEquals(dim(lrr.res[[1]]), c(nscan(blnc), 1))
  sd.res <- medianSdOverAutosomes(lrr.res)
  checkEquals(dim(sd.res), c(nscan(blnc), 2))

  # check error - default incl. values w/o genoData
  checkException(sdByScanChromWindow(blnc, var="LogRRatio"))
  
  close(blnc)
  close(genonc)
}
