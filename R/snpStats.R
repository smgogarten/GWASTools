asSnpMatrix <- function(genoData, snpNames="snpID", scanNames="scanID",
                        snp=c(1,-1), scan=c(1,-1)) {
  
  if (!(requireNamespace("snpStats"))) {
    stop("please install snpStats to use this function")
  }
  
  snpstart <- snp[1]
  snpcount <- snp[2]
  if (snpcount == -1) {
    snpInd <- snpstart:nsnp(genoData)
  } else {
    snpInd <- snpstart:(snpstart+snpcount-1)
  }
  scanstart <- scan[1]
  scancount <- scan[2]
  if (scancount == -1) {
    scanInd <- scanstart:nscan(genoData)
  } else {
    scanInd <- scanstart:(scanstart+scancount-1)
  }
    
  geno <- getGenotype(genoData, snp=snp, scan=scan)
  
  snp.names <- getSnpVariable(genoData, snpNames, index=snpInd)
  scan.names <- getScanVariable(genoData, scanNames, index=scanInd)

  # convert from (NA,0,1,2) to (0,1,2,3)
  geno2 <- geno
  geno2[is.na(geno)] <- 0 # missing
  geno2[geno %in% 2] <- 1 # AA
  geno2[geno %in% 1] <- 2 # AB/BA
  geno2[geno %in% 0] <- 3 # BB
  # SnpMatrix has (scan, snp)
  geno2 <- t(geno2)
  dimnames(geno2) <- list(scan.names, snp.names)
  mode(geno2) <- "raw"
  
  snpmat <- new("SnpMatrix", geno2)
  return(snpmat)
}
