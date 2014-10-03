test_BAFfromClusterMeans <- function() {
  xyfile <- system.file("extdata", "illumina_qxy.nc", package="GWASdata")
  xyNC <- NcdfIntensityReader(xyfile)
  data(illuminaSnpADF)
  snpAnnot <- illuminaSnpADF
  xyData <- IntensityData(xyNC, snpAnnot=snpAnnot)
  
  # fake ncdf file
  blfile <- tempfile()
  BAFfromClusterMeans(xyData, blfile, file.type="gds", verbose=FALSE)

  # read output
  bl <- GdsIntensityReader(blfile)
  baf <- getBAlleleFreq(bl)
  lrr <- getLogRRatio(bl)

  checkIdentical(c(0,1), range(baf, na.rm=TRUE))
  
  close(xyNC)
  close(bl)
  file.remove(blfile)
}
