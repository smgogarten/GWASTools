
ncdfSubsetCheck <- function(
	parent.ncdf,	# name of the parent netCDF file
	sub.ncdf,		# name of the sub netCDF file to check
	sample.include=NULL,	# vector of sampleIDs for samples that should be in sub.ncdf
        snp.include=NULL, # vector of snpIDs for snps that should be in in sub.ncdf
        verbose=TRUE
) {

	# this assumes that "sampleID" is the only 1D="sample" variable in the netCDF
	
	nc <- open.ncdf(parent.ncdf)
	subnc <- open.ncdf(sub.ncdf)

	# check sampleID
	psid <- get.var.ncdf(nc, "sampleID")
	ssid <- get.var.ncdf(subnc, "sampleID")
	chk <- all(is.element(ssid, psid))
	if(!chk) stop("sampleID in sub netCDF is not a subset of parent netCDF")    
	if (is.null(sample.include)) {
                sample.include <- psid
        }
	chk <- all(is.element(sample.include,ssid) & is.element(ssid,sample.include))
	if(!chk) stop("samples in sub netCDF are not the same as sample.include")

        # check snpID
        psnpid <- nc$dim$snp$vals
        ssnpid <- subnc$dim$snp$vals
	chk <- all(is.element(ssnpid, psnpid))
	if(!chk) stop("snpID in sub netCDF is not a subset of parent netCDF")
	if (is.null(snp.include)) {
                snp.include <- psnpid
        }
	chk <- all(is.element(snp.include,ssnpid) & is.element(ssnpid,snp.include))
	if(!chk) stop("snps in sub netCDF are not the same as snp.include")

	# get variable names and their dimensions from subnc
	m <- subnc$nvars
	svardim <- list()
	for(j in 1:m){  # for each variable
		vname <- subnc$var[[j]]$name
		nd <- length(subnc$var[[j]]$dim)   	# number of dimensions
		dname <- rep(NA,nd)
		for(k in 1:nd) dname[k] <- subnc$var[[j]]$dim[[k]]$name		# names of dimensions
		svardim[[j]] <- dname
		names(svardim)[j] <- vname
	}

	# get variable names and their dimensions from parent nc
	n <- nc$nvars
	vardim <- list()
	for(j in 1:n){  # for each variable
		vname <- nc$var[[j]]$name
		nd <- length(nc$var[[j]]$dim)   	# number of dimensions
		dname <- rep(NA,nd)
		for(k in 1:nd) dname[k] <- nc$var[[j]]$dim[[k]]$name		# names of dimensions
		vardim[[j]] <- dname
		names(vardim)[j] <- vname
	}

	# compare number and dimension types of variables
	if(n!=m) stop("number of variables is not the same")
	chk <- list()
	for(i in 1:m) chk[[i]] <- svardim[[i]]==vardim[[i]]
	chk <- all(unlist(chk))
	if(!chk) stop("variables and or their dimension types are not the same")

	# logical vector for selecting samples
	snpsel <- is.element(psnpid, snp.include)

	# check the 1D = "snp" variables
	for(j in 1:m) {
		vname <- names(vardim)[j]
		dname <- vardim[[j]]
		chk <- length(dname)==1 & is.element("snp",dname)
		if(!chk) {
                	next
		} else {
			dat <- get.var.ncdf(nc, vname)[snpsel]
			sdat <- get.var.ncdf(subnc, vname)
			if(length(dat)!=length(sdat)) stop(paste("lengths of variable", vname, "are not the same"))
			chk <- all(is.na(dat)==is.na(sdat))
			if(!chk) stop(paste("elements missing for variable", vname, "are not the same"))
			chk <- all(dat[!is.na(dat)]==sdat[!is.na(dat)])
			if(!chk) stop(paste("values of variable", vname, "are not the same"))
		}
	}

	# check 2D variables
	ns <- length(ssid)
	for(j in 1:m) {
		vname <- names(vardim)[j]
		dname <- vardim[[j]]
		chk <- length(dname)==2 & is.element("sample",dname)
		if(!chk) { 
        		next
		} else {
			for(k in 1:ns){
				if (verbose & k%%10==0) message(paste("checking", vname, "sample", k, "of", ns))
				sdat <- get.var.ncdf(subnc, vname, start=c(1,k), count=c(-1,1))
				index <- which(is.element(psid, ssid[k]))
				pdat <- get.var.ncdf(nc, vname, start=c(1,index), count=c(-1,1))[snpsel]
				chk <- all(is.na(sdat)==is.na(pdat))
				if(!chk) stop(paste("elements missing for variable", vname, "are not the same"))
				chk <- all(sdat[!is.na(sdat)]==pdat[!is.na(sdat)])
				if(!chk) stop(paste("values for variable", vname, "are not the same"))
			}
		}
	}

	message("congratulations, all variables check out")
}



		
