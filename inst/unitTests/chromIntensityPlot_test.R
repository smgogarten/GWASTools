test_chromIntensityPlot <- function() {
  data(illuminaScanADF)
  scanAnnot <- illuminaScanADF
  data(illuminaSnpADF)
  snpAnnot <- illuminaSnpADF
  
  blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
  blnc <- NcdfIntensityReader(blfile)
  blData <-  IntensityData(blnc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

  xyfile <- system.file("extdata", "illumina_qxy.nc", package="GWASdata")
  xync <- NcdfIntensityReader(xyfile)
  xyData <-  IntensityData(xync, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

  genofile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
  genonc <- NcdfGenotypeReader(genofile)
  genoData <-  GenotypeData(genonc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

  scanID <- getScanID(scanAnnot, index=1)
  
  chromIntensityPlot(scan.ids=scanID, chrom.ids=22,
    intenData=blData, type="BAF/LRR", 
    colorGenotypes=TRUE, genoData=genoData)
  
  chromIntensityPlot(scan.ids=scanID, chrom.ids=23,
    intenData=blData, type="BAF", code="interesting sample",
    colorBatch=TRUE, batch.column="BeadSetID")
  
  chromIntensityPlot(scan.ids=scanID, chrom.ids=22,
    intenData=blData, type="LRR", main.txt="interesting sample")
  
  chromIntensityPlot(scan.ids=scanID, chrom.ids=22,
    intenData=xyData, type="R")
  
  chromIntensityPlot(scan.ids=scanID, chrom.ids=22,
    intenData=xyData, type="Theta")
  
  chromIntensityPlot(scan.ids=scanID, chrom.ids=22,
    intenData=xyData, type="R/Theta")

  excl <- getSnpID(blData, index=(getChromosome(blData) == 22 & getPosition(blData) > 3e7 & getPosition(blData) < 4e7))
  chromIntensityPlot(scan.ids=scanID, chrom.ids=22,
    intenData=blData, type="BAF/LRR", snp.exclude=excl)
                         
  checkException(chromIntensityPlot(scan.ids=scanID, chrom.ids=22,
                                        intenData=xyData, type="foo"))
  
  close(genoData)
  close(xyData)
  close(blData)
}
