test_BAFfromGenotypes <- function() {
  data(affySnpADF)
  snpAnnot <- affySnpADF
  data(affyScanADF)
  scanAnnot <- affyScanADF
  xyfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  xyNC <- NcdfIntensityReader(xyfile)
  xyData <- IntensityData(xyNC, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  genofile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
  genoNC <- NcdfGenotypeReader(genofile)
  genoData <- GenotypeData(genoNC, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  
  # fake ncdf file
  blfile <- tempfile()
  BAFfromGenotypes(xyData, genoData, blfile, file.type="gds",
                   call.method="by.plate", plate.name="plate", verbose=FALSE)
  blfile2 <- tempfile()
  BAFfromGenotypes(xyData, genoData, blfile2, file.type="ncdf",
                   call.method="by.plate", plate.name="plate", verbose=FALSE)

  # read output
  bl <- GdsIntensityReader(blfile)
  blData <- IntensityData(bl, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  blNC <- NcdfIntensityReader(blfile2)
  blData2 <- IntensityData(blNC, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  baf <- getBAlleleFreq(blData)
  lrr <- getLogRRatio(blData)
  baf2 <- getBAlleleFreq(blData2)
  lrr2 <- getLogRRatio(blData2)
  checkEquals(baf, baf2)
  checkEquals(lrr, lrr2)

  checkIdentical(c(0,1), range(baf, na.rm=TRUE))
  close(blData)
  close(blData2)

  # by study
  BAFfromGenotypes(xyData, genoData, blfile, file.type="ncdf",
                   call.method="by.study")

  # read output
  blNC <- NcdfIntensityReader(blfile)
  baf <- getBAlleleFreq(blNC)
  lrr <- getLogRRatio(blNC)

  checkIdentical(c(0,1), range(baf, na.rm=TRUE))
  
  close(xyData)
  close(genoData)
  close(blNC)
  file.remove(blfile, blfile2)
}
