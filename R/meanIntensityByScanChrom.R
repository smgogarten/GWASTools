
meanIntensityByScanChrom <- function(intenData, # object of type IntensityData
                                           vars = c("X", "Y"),
                                           snp.exclude = NULL,
                                           verbose = TRUE)
{
        # check number of variables
        nVars <- length(vars)
        if (nVars > 2) stop("please specify one or two variables")

        # check that variables are present in data
        for (v in vars) {
                if (!hasVariable(intenData, v)) stop(paste(v, "not found in intenData"))
        }
        
	# get start and count row indices for each chromosome type
        chrom <- getChromosome(intenData, char=TRUE)
        uniqChrom <- unique(chrom)
        nChrom <- length(uniqChrom)
	indices <- matrix(NA, nChrom, 2, dimnames=list(uniqChrom, c("start","stop")))
	for(i in 1:nChrom) indices[i,] <- range(which(is.element(chrom, uniqChrom[i])))

	# logical vectors for excluding snps by chromosome
        snpID <- getSnpID(intenData)
        nSnp <- length(snpID)
	if(length(snp.exclude)!=0) {
		exclude <- !is.element(snpID, snp.exclude)
		chrex <- vector("list", nChrom)
		for(j in 1:nChrom){ 
			chrex[[j]] <- exclude[indices[j,"start"]:indices[j,"stop"]]
		}
	}

        scanID <- getScanID(intenData)
        nScan <- length(scanID)
        
	# matrices to hold mean and sd of intensities or other quantity variable
	mnx <- matrix(NA, nScan, nChrom, dimnames=list(scanID, uniqChrom))
	sdx <- matrix(NA, nScan, nChrom, dimnames=list(scanID, uniqChrom))
	if (nVars > 1) {
                mny <- matrix(NA, nScan, nChrom, dimnames=list(scanID, uniqChrom)) 
                sdy <- matrix(NA, nScan, nChrom, dimnames=list(scanID, uniqChrom))
                mni <- matrix(NA, nScan, nChrom, dimnames=list(scanID, uniqChrom))
                sdi <- matrix(NA, nScan, nChrom, dimnames=list(scanID, uniqChrom))
        }
                
	for(i in 1:nScan){	# for each sample
		if (verbose & (i %% 100 == 0))
			message(paste("scan",i,"of",nScan))
                
                x <- getVariable(intenData, vars[1], snp=c(1,-1), scan=c(i,1))
                        
                if (nVars > 1) {
                         y <- getVariable(intenData, vars[2], snp=c(1,-1), scan=c(i,1))
                         inten <- x + y
                }
                
                for (j in 1:nChrom) {	# for each chromosome
                        xt <- x[indices[j,"start"]:indices[j,"stop"]]
                        # remove snps to be excluded
                        if(length(snp.exclude)!=0) xt <- xt[chrex[[j]]]
                                
                        if (nVars > 1) {
                                yt <- y[indices[j,"start"]:indices[j,"stop"]]
                                if(length(snp.exclude)!=0) yt <- yt[chrex[[j]]]

                                # for 2 variables, remove NAs in either variable
                                xd <- xt[!is.na(xt) & !is.na(yt)]
                                yd <- yt[!is.na(xt) & !is.na(yt)]

                                # mean and sd for y
                                mny[i,j] <- mean(yd, na.rm=TRUE)   
                                sdy[i,j] <- sd(yd, na.rm=TRUE)                                 
                                
                                # intensity for the jth chromosome
                                z <- inten[indices[j,"start"]:indices[j,"stop"]]
                                if(length(snp.exclude)!=0) z <- z[chrex[[j]]]
                                mni[i,j] <- mean(z, na.rm=TRUE)
                                sdi[i,j] <- sd(z, na.rm=TRUE)
                        } else {
                                # nothing to change
                                xd <- xt
                        }

                        # mean and SD for x
                        mnx[i,j] <- mean(xd, na.rm=TRUE)
                        sdx[i,j] <- sd(xd, na.rm=TRUE)      
                                
                        rm(xt,xd)
                        if (nVars > 1) rm(z,yt,yd)
                }
                rm(x)
                if (nVars > 1) rm(y,inten)
        }

        if (nVars > 1) {
                r <- list(mni,sdi,mnx,sdx,mny,sdy)
                names(r) <- c("mean.intensity", "sd.intensity",
                              paste(c("mean","sd"), vars[1], sep="."),
                              paste(c("mean","sd"), vars[2], sep="."))
        } else {
                r <- list(mnx,sdx)
                names(r) <- paste(c("mean","sd"), vars[1], sep=".")
        }
	return(r)
}
