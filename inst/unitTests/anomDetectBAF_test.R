check_anomSegmentBAF <- function(blData, genoData) {
  scan.ids <- getScanID(blData, index=1:2)
  chrom.ids <- unique(getChromosome(blData))
  snp.ids <- getSnpID(blData)
  seg <- anomSegmentBAF(blData, genoData, scan.ids=scan.ids,
                        chrom.ids=chrom.ids, snp.ids=snp.ids)
  checkTrue(all(seg$scanID %in% scan.ids))
  checkTrue(all(seg$chrom %in% chrom.ids))
  return(seg)
}

check_anomFilterBAF <- function(blData, genoData, segments) {
  snp.ids <- getSnpID(blData)
  data(centromeres.hg18)
  filt <- anomFilterBAF(blData, genoData, segments=segments, snp.ids=snp.ids,
                        centromere=centromeres.hg18)
  checkEquals(segments, filt$raw[,1:6])
  checkTrue(all(filt$filtered$scanID %in% segments$scanID))
  checkTrue(all(filt$filtered$chrom %in% segments$chrom))
  return(filt)
}

check_anomDetectBAF <- function(blData, genoData) {
  scan.ids <- getScanID(blData, index=1:2)
  chrom.ids <- unique(getChromosome(blData))
  snp.ids <- getSnpID(blData)
  data(centromeres.hg18)
  anom <- anomDetectBAF(blData, genoData, scan.ids=scan.ids,chrom.ids=chrom.ids,
                        snp.ids=snp.ids, centromere=centromeres.hg18,
                        num.mark.thresh=20)
  anom <- anomDetectBAF(blData, genoData, scan.ids=scan.ids,chrom.ids=chrom.ids,
                        snp.ids=snp.ids, centromere=centromeres.hg18)
  checkEquals(length(anom), 4)
  return(anom)
}

test_anomDetectBAF <- function() {
  data(illuminaScanADF)
  scanAnnot <- illuminaScanADF
  data(illuminaSnpADF)
  snpAnnot <- illuminaSnpADF
  
  blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
  blnc <- NcdfIntensityReader(blfile)
  blData <-  IntensityData(blnc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

  genofile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
  genonc <- NcdfGenotypeReader(genofile)
  genoData <-  GenotypeData(genonc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

  seg <- check_anomSegmentBAF(blData, genoData)
  filt <- check_anomFilterBAF(blData, genoData, seg)
  anom <- check_anomDetectBAF(blData, genoData)
  checkEquals(filt, anom)
  
  close(blData)
  close(genoData)
}
