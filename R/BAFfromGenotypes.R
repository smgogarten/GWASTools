
## since we calculate BAF/LRR by SNP, need to add data by snp
.createGdsBySnp <- function(sample.id, snp.annotation, filename, variables, precision,
                       compress) {

    ## define precision for gds
    precision <- ifelse(precision == "double", "float64", "float32")

    ## create gds file
    gds <- createfn.gds(filename)

    ## add standard variables
    add.gdsn(gds, "sample.id", sample.id, compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.id", snp.annotation$snpID, compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.chromosome", snp.annotation$chromosome, storage="uint8",
             compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.position", snp.annotation$position, compress=compress, closezip=TRUE)
    sync.gds(gds)

    ## add selected variables
    for (v in variables) {
        add.gdsn(gds, v, storage=precision, valdim=c(nrow(snp.annotation), length(sample.id)))
    }

    sync.gds(gds)
    gds
}

.addDataBySnp <- function(x, ...) UseMethod(".addDataBySnp", x)
.addDataBySnp.gds.class <- function(x, vars, dat, snp.start, snp.count) {
    for (v in vars) {
        write.gdsn(index.gdsn(x, v), val=dat[[v]], start=c(snp.start,1), count=c(snp.count,-1))
    }
}

.addDataBySnp.ncdf <- function(x, vars, dat, snp.start, snp.count) {
    for (v in vars) {
        put.var.ncdf(x, v, vals=dat[[v]], start=c(snp.start,1), count=c(snp.count,-1))
    }
}

BAFfromGenotypes <- function(
		intenData, 
		genoData, 
		filename,
                file.type=c("gds", "ncdf"),
		min.n.genotypes = 2, 
		call.method = c("by.plate", "by.study"),
		plate.name = "plate",
		block.size = 5000, 
                precision = "single", compress = "ZIP_RA",
		verbose = TRUE) {
				
  # check that dimensions of intenData and genoData are equal
  intenSnpID <- getSnpID(intenData)
  genoSnpID <- getSnpID(genoData)
  if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
  intenScanID <- getScanID(intenData)
  genoScanID <- getScanID(genoData)
  if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")

  # Check that plate is included in data if by.plate is being used
  call.method <- match.arg(call.method)
  if(call.method == "by.plate") {
    if (hasScanVariable(intenData, plate.name)) {
      plate <- getScanVariable(intenData, plate.name)
    } else if (hasScanVariable(genoData, plate.name)) {
      plate <- getScanVariable(genoData, plate.name)
    } else stop("call.method==by.plate but plate.name not found in intenData or genoData")
  }
	
  ## get file type
  file.type <- match.arg(file.type)

  ## create data file
  snp.annotation <- getSnpVariable(intenData, c("snpID", "chromosome", "position"))
  variables <- c("BAlleleFreq", "LogRRatio")
  if (file.type == "gds") {
      genofile <- .createGdsBySnp(intenScanID, snp.annotation, filename, variables, precision, compress)
  } else if (file.type == "ncdf") {
      genofile <- .createNcdf(snp.annotation, filename, variables, nscan(intenData),
                              precision)
      put.var.ncdf(genofile, "sampleID", vals=intenScanID)
  }

  # make a matrix of T/F indicating which plates the samples are on
  N <- length(intenScanID)
  if(call.method == "by.plate") {
    sidPlate <- cbind(intenScanID, plate)
    orderedPlate <- levels(as.factor(plate))
    Q <- length(orderedPlate)

    P <- matrix(rep(FALSE,Q*N), nrow=Q, ncol=N)
    dimnames(P) <- list(orderedPlate, intenScanID)
    for(i in 1:nrow(sidPlate)) {
      P[sidPlate[i,2],sidPlate[i,1]] <- TRUE
    }
  } else {
    # fill P with all trues to essentially make a one plate study if want to analyze 'by study' rather than 'by plate'
    Q <- 1
    P <- matrix(rep(TRUE, Q*N), nrow=Q, ncol=N)
  }

  # define block.sizes
  n <- length(intenSnpID)
  nblock <- floor(n/block.size)
  remain <- n-block.size*nblock
  rm(n)
  if(remain == 0) {M <- nblock}  else { M <- nblock+1 }
  # M is the number of iterations needed to get all samples

  m <- 1 # first snp of first block.size

  # loop through block.sizes, get X, Y and G and use them to calculate B and L

  for (i in 1:M) {
	if (verbose)
		message(paste("block", i, "of", M))
	if(i <= nblock) { n <- block.size } else { n <- remain }
	# size of current block.size = number of snps to get		# CHANGE columns to snps (rows)

        geno <- getGenotype(genoData, snp=c(m,n), scan=c(1,-1), drop=FALSE)
        x <- getX(intenData, snp=c(m,n), scan=c(1,-1), drop=FALSE)
        y <- getY(intenData, snp=c(m,n), scan=c(1,-1), drop=FALSE)

        # output matrices
        BAF <- matrix(0.0, nrow=n, ncol=N)
        LRR <- matrix(0.0, nrow=n, ncol=N)
        
	for (s in 1:n) { # for each snp within the block.size

		# get geno, x and y vectors for the current snp
		gv <- geno[s,]; xv <- x[s,]; yv <- y[s,]

		T <- rep(NA, N); R <- rep(NA, N)
		B <- rep(NA, N); L <- rep(NA, N)	
		
		T <- atan(yv/xv)*(2/pi)		
		R <- xv+yv

		for(k in 1:Q) { # loop through plates, if by study, essentially one plate so Q=1
		
			p <- P[k,] # T for all sample.ids that are on the current plate k
		#	sum(p) # check to ensure this is the number of samples per plate
			
			AAind <- p & gv==2 & !is.na(gv)	
			ABind <- p & gv==1 & !is.na(gv)	
			BBind <- p & gv==0 & !is.na(gv)			
			
			if(sum(AAind) < min.n.genotypes | sum(ABind) < min.n.genotypes | sum(BBind) < min.n.genotypes)
			{ next }

			tAA <- mean(T[AAind], na.rm=TRUE)		# ADDED na.rm=TRUE arg
			tAB <- mean(T[ABind], na.rm=TRUE)		# note:  can't use na.rm=T because T is theta
			tBB <- mean(T[BBind], na.rm=TRUE)
							
			rAA <- mean(R[AAind], na.rm=TRUE)
			rAB <- mean(R[ABind], na.rm=TRUE)
			rBB <- mean(R[BBind], na.rm=TRUE)

			TC <- matrix(rep(NA, 4*N),nrow=4, ncol=N)
			TC[1,] <- T < tAA
			TC[2,] <- tAA <= T & T < tAB
			TC[3,] <- tAB <= T & T < tBB
			TC[4,] <- tBB <= T

			# check to ensure there are no NAs in TC matrix
			# TC[is.na(TC)] # should be no entries			# WHAT IF THERE ARE NAs?
			TC[is.na(TC)] <- F						# ADDED

			# calculate and store B allele frequencies
			# ensure only loading for those samples on the current plate
			B[p & TC[1,]] <- 0
			c2 <- p & TC[2,]
			B[c2] <- (.5*(T[c2]-tAA))/(tAB-tAA)
			c3 <- p & TC[3,]
			B[c3] <- .5+(.5*(T[c3]-tAB))/(tBB-tAB)
			B[p & TC[4,]] <- 1

			# calculate log R ratio for samples on current plate
			c <- (p & TC[1,]) | (p & TC[2,])
			rhat <- rAA + ((T[c]-tAA)*(rAB-rAA))/(tAB-tAA)	# if rhat is negative, make it missing (NA)
			rhat[which(rhat<0)] <- NA
			L[c] <- log2(R[c]/rhat)

			c2 <- (p & TC[3,]) | (p & TC[4,])
			rhat2 <- rAB + ((T[c2]-tAB)*(rBB-rAB))/(tBB-tAB)	# if rhat2 is negative, make it missing (NA)
			rhat2[which(rhat2<0)] <- NA
			L[c2] <- log2(R[c2]/rhat2)

		} # end for loop through plates
                BAF[s,] <- B
                LRR[s,] <- L
		
	} # end for loop through SNPs within block.size

        dat <- list("BAlleleFreq"=BAF, "LogRRatio"=LRR)
        .addDataBySnp(genofile, variables, dat, snp.start=m, snp.count=n)

	m <- m+n

      } # end for loop through number of block.sizes

  .close(genofile, verbose=verbose)
  return(invisible(NULL))
}
			


	


