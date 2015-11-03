
BAFfromClusterMeans <- function(intenData,
                                filename, file.type=c("gds", "ncdf"),
                                clusterMeanVars = c("tAA","tAB","tBB","rAA","rAB","rBB"),
                                precision = "single", compress = "ZIP_RA",
                                verbose = TRUE)
{
  # check that cluster means are in intenData
  stopifnot(all(hasSnpVariable(intenData, clusterMeanVars)))

  # get cluster means
  tAA <- getSnpVariable(intenData, clusterMeanVars[1])
  tAB <- getSnpVariable(intenData, clusterMeanVars[2])
  tBB <- getSnpVariable(intenData, clusterMeanVars[3])
  rAA <- getSnpVariable(intenData, clusterMeanVars[4])
  rAB <- getSnpVariable(intenData, clusterMeanVars[5])
  rBB <- getSnpVariable(intenData, clusterMeanVars[6])
  
  # find line segments connecting AA,AB and AB,BB
  slopeA <- (rAB - rAA) / (tAB - tAA)
  interceptA <- rAB - slopeA*tAB
  slopeB <- (rBB - rAB) / (tBB - tAB)
  interceptB <- rAB - slopeB*tAB
  
  ## get file type
  file.type <- match.arg(file.type)

  ## create data file
  snp.annotation <- getSnpVariable(intenData, c("snpID", "chromosome", "position"))
  variables <- c("BAlleleFreq", "LogRRatio")
  if (file.type == "gds") {
      genofile <- .createGds(snp.annotation, filename, variables, precision, compress)
  } else if (file.type == "ncdf") {
      genofile <- .createNcdf(snp.annotation, filename, variables, nscan(intenData),
                              precision)
  }


  # loop over samples
  scanID <- getScanID(intenData)
  nScan <- length(scanID)
  nSnp <- length(tAA)
  for (i in 1:nScan) { 
    if (verbose & (i %% 100 == 0)) 
      message(paste("scan", i, "of", nScan))
    
    x <- getX(intenData, snp=c(1,-1), scan=c(i,1))
    y <- getY(intenData, snp=c(1,-1), scan=c(i,1))

    theta <- atan(y/x)*(2/pi)
    r <- x+y

    # divide data into quadrants
    tsec <- matrix(nrow=nSnp, ncol=4)
    tsec[,1] <- theta < tAA
    tsec[,2] <- tAA <= theta & theta < tAB
    tsec[,3] <- tAB <= theta & theta < tBB
    tsec[,4] <- tBB <= theta
    tsec[is.na(tsec)] <- FALSE

    # calculate theta distance along line segment
    tfracA <- (theta - tAA) / (tAB - tAA)
    tfracB <- (theta - tAB) / (tBB - tAB)

    # calculate BAF
    baf <- rep(NA, nSnp)
    baf[tsec[,1]] <- 0
    baf[tsec[,2]] <- 0.5*tfracA[tsec[,2]]
    baf[tsec[,3]] <- 0.5*tfracB[tsec[,3]] + 0.5
    baf[tsec[,4]] <- 1
    
    # calculate R_expected
    rexp <- rep(NA, nSnp)
    rexpA <- slopeA*theta + interceptA
    rexp[tsec[,1] | tsec[,2]] <- rexpA[tsec[,1] | tsec[,2]]
    rexpB <- slopeB*theta + interceptB
    rexp[tsec[,3] | tsec[,4]] <- rexpB[tsec[,3] | tsec[,4]]
    # if rexp is negative, make it NA
    rexp[rexp < 0] <- NA
    
    # log R Ratio
    lrr <- log2(r/rexp)

    # write to file
    dat <- list("BAlleleFreq"=baf, "LogRRatio"=lrr)
    .addData(genofile, variables, dat, scanID[i], i)
  }
  .close(genofile, verbose=verbose)
  return(invisible(NULL))
}
