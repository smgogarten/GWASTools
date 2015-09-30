genoClusterPlotByBatch <- function(intenData, 
                         genoData, 
                         plot.type=c("RTheta","XY"),
                         snpID, 
                         batchVar,
                         main.txt = NULL,
                         scan.sel = NULL, 
                         colors = c("default", "neon", "primary"),
                         verbose = TRUE,
                         ...)
{

  if (verbose) message(paste("start time =", Sys.time()))
  colors <- match.arg(colors)
			
  # check that dimensions of intenData and genoData are equal
  intenSnpID <- getSnpID(intenData)
  genoSnpID <- getSnpID(genoData)
  if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
  intenScanID <- getScanID(intenData)
  genoScanID <- getScanID(genoData)
  if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
    
  # Get snp.index corresponding to each snpID (while keeping same order as snpID and main.txt)
  if(!all(is.element(snpID, intenSnpID))) {
    stop("At least one selected SNP in snpID is invalid, make sure to use integer snpID")
  }
  nsnp <- length(snpID)
  snp.index <- rep(NA,nsnp)
  for(i in 1:nsnp) { 
    snp.index[i] <- which(is.element(intenSnpID, snpID[i])) 
  }	

  # get plot type
  plot.type <- match.arg(plot.type)
  
  # Get sample index corresponding to each scan.sel to include in plot
  if(!is.null(scan.sel)) {
    if(!all(is.element(scan.sel, intenScanID))) {
      stop("At least one selected sample in scan.sel is not in the genotype file")
    }
    samp.plot <- is.element(intenScanID, scan.sel)
  } else {
    samp.plot <- rep(TRUE, length(intenScanID))
  }

  # check that batchVar is present
  if (hasScanVariable(intenData, batchVar)) {
    sample.ann.plates <- getScanVariable(intenData, batchVar)
  } else  if (hasScanVariable(genoData, batchVar)) {
    sample.ann.plates <- getScanVariable(genoData, batchVar)
  } else stop ("batchVar not found in intenData or genoData")
  plates <- unique(sample.ann.plates[samp.plot])
  np <- length(plates);

  if (is.null(main.txt)) {
    main.txt <- paste("snp", as.character(snpID))
  }

  for(i in 1:nsnp) { # loop through snp int.ids
    geno <- getGenotype(genoData, snp=c(snp.index[i],1), scan=c(1,-1))	
    y <- getY(intenData, snp=c(snp.index[i],1), scan=c(1,-1))
    x <- getX(intenData, snp=c(snp.index[i],1), scan=c(1,-1))
    for(j in 1:np){ # loop through plates
      cplate <- plates[j];
      samp.logic1 <- sample.ann.plates == cplate
      if(!is.null(scan.sel)) {
        samp.logic <- samp.logic1 & samp.plot
      } else { samp.logic <- samp.logic1 }
      if (verbose) message(paste("i=",i,"j=",j, Sys.time()))
      genosm <- geno[samp.logic];
      ysm <- y[samp.logic];
      xsm <- x[samp.logic];
      stopifnot(length(genosm)==length(ysm) & length(ysm)==length(xsm))
      xcol <- rep(NA, length(genosm))
      cols <- .colorByGeno(colors)
      xcol[is.na(genosm)] <- cols["NA"]
      xcol[genosm==0] <- cols["BB"]
      xcol[genosm==1] <- cols["AB"]
      xcol[genosm==2] <- cols["AA"]
      stopifnot(sum(is.na(xcol))==0)
      xpch <- rep(NA, length(genosm))
      xpch[is.na(genosm)] <- 4
      xpch[!is.na(genosm)] <- 1
      stopifnot(sum(is.na(xpch))==0)
      txt <- paste(main.txt[i], "\n", as.character(cplate), sep="")
      if(plot.type=="RTheta") {
        theta <- atan(ysm/xsm)*(2/pi)
        r <- xsm+ysm
        plot(theta, r, xlab="Theta", ylab="R", xlim=c(0,1),col=xcol, pch=xpch, main=txt, cex=0.8, ...)
      } else { 
        plot(xsm, ysm, xlab="X", ylab="Y",col=xcol, pch=xpch, main=txt, cex=0.8, ...)
      }
      rm(genosm); rm(xsm); rm(ysm)
    } #j
  }  #i
}
