genoClusterPlot <- function(intenData, 
                         genoData, 
                         plot.type=c("RTheta","XY"),
                         snpID, 
                         main.txt = NULL,
                         by.sex= FALSE, # if by.sex is TRUE, sex annotation should be included in intenData or genoData
                         scan.sel = NULL, 
                         scan.hilite = NULL,
                         start.axis.at.0 = FALSE,
                         verbose = TRUE,
                         ...)
{

  if (verbose) message(paste("start time =", Sys.time()))
			
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
  nplot <- length(snpID)
  snp.index <- rep(NA,nplot)
  for(i in 1:nplot) { 
    snp.index[i] <- which(is.element(intenSnpID, snpID[i])) 
  }	

  # get plot type
  plot.type <- match.arg(plot.type)
  
  # Get sample index corresponding to each scan.sel to include in plot
  if(!is.null(scan.sel)) {
    if(!all(is.element(scan.sel, intenScanID))) {
      stop("At least one selected sample in scan.sel is not in the genotype file")
    }
    scan.index <- which(is.element(intenScanID, scan.sel))
  } else {
    scan.index <- 1:length(intenScanID)
  }
		
  # Get indicator for samples to be hilited in the plots
  if(!is.null(scan.hilite)) {
    if(!all(is.element(scan.hilite, intenScanID))) {
      stop("At least one sample in scan.hilite is not in the genotype file")
    }
    hilite.ind <- ifelse(is.element(intenScanID, scan.hilite),1,0)
    hilite.ind <- hilite.ind[scan.index]
  } else { 
    hilite.ind <- rep(0,length(scan.index))
  }
  	
  # Check that sex is included in data if by.sex is being used
  if(by.sex) {
    if (hasSex(intenData)) {
      sex <- getSex(intenData)
    } else if (hasSex(genoData)) {
      sex <- getSex(genoData)
    } else stop("by.sex=TRUE but sex annotation not found in intenData or genoData")
    by.sex <- sex[scan.index]
  }

  if (is.null(main.txt)) {
    main.txt <- paste("snp", as.character(snpID))
  }
	
  for(i in 1:nplot) {
    if (verbose) 
      message(paste("plot number =",i, " system time =", Sys.time()))

    geno <- getGenotype(genoData, snp=c(snp.index[i],1), scan=c(1,-1))
    geno <- geno[scan.index]	

    y <- getY(intenData, snp=c(snp.index[i],1), scan=c(1,-1))
    y <- y[scan.index] 
    x <- getX(intenData, snp=c(snp.index[i],1), scan=c(1,-1))
    x <- x[scan.index] 
    xcol <- rep(NA, length(geno))   
    xcol[is.na(geno)] <- "black"
    xcol[geno==0 & hilite.ind==0] <- "blue"
    xcol[geno==1 & hilite.ind==0] <- "green"
    xcol[geno==2 & hilite.ind==0] <- "red"
    xcol[geno==0 & hilite.ind==1] <- "orange"
    xcol[geno==1 & hilite.ind==1] <- "magenta"
    xcol[geno==2 & hilite.ind==1] <- "yellow"
    xpch <- rep(NA, length(geno)) 
    xpch[is.na(geno)] <- 4
    xpch[!is.na(geno)] <- 1
    xpch[hilite.ind==1] <- 9
    
    if(sum(hilite.ind)==0) {
      cex <- rep(1, length(geno))
    } else{
      cex <- rep(0.8, length(geno))
      cex[hilite.ind==1] <- 1.5
    }
    # Note: by.sex takes precedence for plotting character if it is used
    # but the colors will still be different for hilited samples
    if(!is.null(by.sex)) {
      xpch[by.sex=="F"] <- 1
      xpch[by.sex=="M"] <- 3
    }

    if(plot.type=="RTheta") {
      theta <- atan(y/x)*(2/pi)
      r <- x+y
      if (start.axis.at.0) {
        plot(theta, r, xlab="Theta", ylab="R", xlim=c(0,1), ylim=c(0,max(r,na.rm=TRUE)), col=xcol, pch=xpch, main=main.txt[i], cex=cex, ...)
      } else {
        plot(theta, r, xlab="Theta", ylab="R", xlim=c(0,1), col=xcol, pch=xpch, main=main.txt[i], cex=cex, ...)
      }
      points(theta[hilite.ind==1],r[hilite.ind==1],col=xcol[hilite.ind==1],pch=xpch[hilite.ind==1], cex=cex[hilite.ind==1], ...)
    } 
    else {
      if (start.axis.at.0) {
        plot(x, y, xlab="X", ylab="Y", xlim=c(0,max(x,na.rm=TRUE)), ylim=c(0,max(y,na.rm=TRUE))), col=xcol, pch=xpch, main=main.txt[i], cex=cex, ...)
      } else { 
        plot(x, y, xlab="X", ylab="Y", col=xcol, pch=xpch, main=main.txt[i], cex=cex, ...)
      }
      points(x[hilite.ind==1],y[hilite.ind==1],col=xcol[hilite.ind==1],pch=xpch[hilite.ind==1], cex=cex[hilite.ind==1], ...)
    }
  }
}
