test_pseudoautoIntensityPlot <- function() {
  data(illumina_scan_annot)
  scanAnnot <- ScanAnnotationDataFrame(illumina_scan_annot)
  
  blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
  blnc <- NcdfIntensityReader(blfile)
  blData <-  IntensityData(blnc, scanAnnot=scanAnnot)

  scanID <- getScanID(scanAnnot, index=1)
  pseudoautoIntensityPlot(scan.ids=scanID, intenData=blData)

  pseudoautoIntensityPlot(scan.ids=scanID, intenData=blData, plotY=TRUE)
  close(blData)
}
