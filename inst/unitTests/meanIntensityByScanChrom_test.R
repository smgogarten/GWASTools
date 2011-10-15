test_meanIntensityByScanChrom <- function() { 
  # simulate data
  ncfile <- tempfile()
  simulateIntensityMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfIntensityReader(ncfile)
  intenData <- IntensityData(nc)
  scanID <- getScanID(intenData)
  chrom <- getChromosome(intenData, char=TRUE)
  uchr <- unique(chrom)
  x <- getX(intenData)
  y <- getY(intenData)

  # expected results
  nchr <- 26
  mn.X <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  mn.Y <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  sd.X <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  sd.Y <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  mn.XY <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  sd.XY <- matrix(NA, nrow=length(scanID), ncol=length(uchr),
                 dimnames=list(scanID, uchr))
  for (i in 1:nchr) {
    xc <- x[(chrom == uchr[i]),]
    yc <- y[(chrom == uchr[i]),]
    mn.X[,i] <- colMeans(xc, na.rm=TRUE)
    mn.Y[,i] <- colMeans(yc, na.rm=TRUE)
    sd.X[,i] <- apply(xc, 2, sd)
    sd.Y[,i] <- apply(yc, 2, function(f) sd(f, na.rm=TRUE))
    mn.XY[,i] <- colMeans(xc+yc, na.rm=TRUE)
    sd.XY[,i] <- apply(xc+yc, 2, function(f) sd(f, na.rm=TRUE))
  }
  exp <- list("mean.intensity"=mn.XY, "sd.intensity"=sd.XY,
              "mean.X"=mn.X, "sd.X"=sd.X, "mean.Y"=mn.Y, "sd.Y"=sd.Y)

  res <- meanIntensityByScanChrom(intenData)
  checkIdentical(exp, res)

  # one var only
  exp <- list("mean.X"=mn.X, "sd.X"=sd.X)
  res <- meanIntensityByScanChrom(intenData, vars="X")
  checkIdentical(exp, res)

  # snp.exclude
  snpID <- getSnpID(intenData)
  snp.exclude <- snpID[c(1,2,11)]
  for (i in 1:nchr) {
    xc <- x[(chrom == uchr[i] & !(snpID %in% snp.exclude)),]
    yc <- y[(chrom == uchr[i] & !(snpID %in% snp.exclude)),]
    mn.X[,i] <- colMeans(xc, na.rm=TRUE)
    mn.Y[,i] <- colMeans(yc, na.rm=TRUE)
    sd.X[,i] <- apply(xc, 2, sd)
    sd.Y[,i] <- apply(yc, 2, function(f) sd(f, na.rm=TRUE))
    mn.XY[,i] <- colMeans(xc+yc, na.rm=TRUE)
    sd.XY[,i] <- apply(xc+yc, 2, function(f) sd(f, na.rm=TRUE))
  }
  exp <- list("mean.intensity"=mn.XY, "sd.intensity"=sd.XY,
              "mean.X"=mn.X, "sd.X"=sd.X, "mean.Y"=mn.Y, "sd.Y"=sd.Y)

  res <- meanIntensityByScanChrom(intenData, snp.exclude=snp.exclude)
  checkIdentical(exp, res)

  close(intenData)                    
  file.remove(ncfile)
}
