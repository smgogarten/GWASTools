missingGenotypeBySnpSex <- function(
		genoData, # object of type GenotypeData
		scan.exclude = NULL, # vector of scanID
		verbose = TRUE) {

  # check that sex in included in genoData
  if (!hasSex(genoData)) stop("sex not found in genoData")

  # select males and females which are not excluded
  scanID <- getScanID(genoData)
  sex <- getSex(genoData)	
  males <- which(!is.na(sex) & sex == "M" & !is.element(scanID, scan.exclude))
  females <- which(!is.na(sex) & sex == "F" & !is.element(scanID, scan.exclude))

  spx <- c(NA,NA)
  names(spx) <- c("M", "F")
  spx["F"] <- length(females)
  spx["M"] <- length(males)

  # matrix to hold missing counts
  snpID <- getSnpID(genoData)
  miss.cnt <- matrix(0, length(snpID), 2, dimnames=list(snpID, c("M","F")))

  N <- length(scanID)
  for (i in 1:N) {
    if (verbose & (i %% 100 == 0)) {
      message(paste("reading scan", i, "of", N))
    }
    if (i %in% males) {
      geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))
      miss.cnt[,"M"] <- miss.cnt[,"M"] + is.na(geno)
    } else if (i %in% females) {
      geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))
      miss.cnt[,"F"] <- miss.cnt[,"F"] + is.na(geno)
    }
  }

  # compute missing fraction
  ychr <- getChromosome(genoData, char=TRUE) == "Y"
  miss.frac <- rowSums(miss.cnt) / sum(spx)
  # for Y chromosome, use only males
  miss.frac[ychr] <- miss.cnt[ychr, "M"] / spx["M"]
  
  miss <- list(miss.cnt, spx, miss.frac)
  names(miss) <- c("missing.counts", "scans.per.sex", "missing.fraction")
  return(miss)			
}
