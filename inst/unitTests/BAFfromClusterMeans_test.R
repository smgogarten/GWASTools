test_BAFfromClusterMeans <- function() {
  xyfile <- system.file("extdata", "illumina_qxy.nc", package="GWASdata")
  xyNC <- NcdfIntensityReader(xyfile)
  data(illuminaSnpADF)
  snpAnnot <- illuminaSnpADF
  xyData <- IntensityData(xyNC, snpAnnot=snpAnnot)
  nsamp <- length(getScanID(xyData))
  
  # fake ncdf file
  blfile <- tempfile()
  ncdfCreate(pData(snpAnnot), blfile, variables=c("BAlleleFreq","LogRRatio"),
                  n.samples=nsamp)

  BAFfromClusterMeans(xyData, blfile, verbose=FALSE)

  # read output
  blNC <- NcdfIntensityReader(blfile)
  blData <- IntensityData(blNC)
  baf <- getBAlleleFreq(blData)
  lrr <- getLogRRatio(blData)

  checkIdentical(c(0,1), range(baf, na.rm=TRUE))
  
  close(xyNC)
  close(blNC)
  file.remove(blfile)
}
