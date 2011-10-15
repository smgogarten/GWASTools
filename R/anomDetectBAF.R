# wrapper function to segment and filter anomalies in one step

anomDetectBAF <- function(intenData, genoData, scan.ids, chrom.ids, snp.ids,
                          centromere, low.qual.ids=NULL, ...) {
  segments <- anomSegmentBAF(intenData, genoData, scan.ids, chrom.ids, snp.ids)
  anoms <- anomFilterBAF(intenData, genoData, segments, snp.ids,
                         centromere, low.qual.ids, ...)
  return(anoms)
}
