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
  nsamp <- length(getScanID(genoData))
  
  # fake ncdf file
  blfile <- tempfile()
  ncdfCreate(affy_snp_annot, blfile, variables=c("BAlleleFreq","LogRRatio"),
                  n.samples=nsamp)

  BAFfromGenotypes(xyData, genoData, blfile,
                   call.method="by.plate", plate.name="plate")

  # read output
  blNC <- NcdfIntensityReader(blfile)
  blData <- IntensityData(blNC, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  baf <- getBAlleleFreq(blData)
  lrr <- getLogRRatio(blData)

  checkIdentical(c(0,1), range(baf, na.rm=TRUE))
  close(blData)

  # by study
  BAFfromGenotypes(xyData, genoData, blfile,
                   call.method="by.study")

  # read output
  blNC <- NcdfIntensityReader(blfile)
  blData <- IntensityData(blNC)
  baf <- getBAlleleFreq(blData)
  lrr <- getLogRRatio(blData)

  checkIdentical(c(0,1), range(baf, na.rm=TRUE))
  
  close(xyData)
  close(genoData)
  close(blData)
  file.remove(blfile)
}
