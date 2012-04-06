
# This function replicates an existing netCDF, but includes only a subset of samples and snps

ncdfSubset <- function(
	parent.ncdf,	# name of the parent netCDF file
	sub.ncdf,		# name of the sub netCDF file to create
	sample.include=NULL,	# vector of sampleIDs for samples to include in sub.ncdf
        snp.include=NULL, # vector of snpIDs for snps to include in sub.ncdf
        verbose=TRUE
) {

	# This function only works for netCDF files having up to two dimensions, named "snp" and "sample"
	nc <- open.ncdf(parent.ncdf)   
	
	# check that sample.include are all elements of "sampleID"
	sampid <- get.var.ncdf(nc,"sampleID")     
	if (is.null(sample.include)) {
                sample.include <- sampid
        }
	chk <- all(is.element(sample.include,sampid))
	if(!chk) stop("sample.include elements are not all members of netCDF sampleID")

	# logical vector for selecting samples
	sampsel <- is.element(sampid, sample.include)

	# check that snp.include are all elements of "snp"
	snpid <- nc$dim$snp$vals   
	if (is.null(snp.include)) {
                snp.include <- snpid
        }
	chk <- all(is.element(snp.include,snpid))
	if(!chk) stop("snp.include elements are not all members of netCDF snp")

	# logical vector for selecting samples
	snpsel <- is.element(snpid, snp.include)

	# Define dimensions
	n <- nc$ndims
	dims <- vector("list",n)
	for(i in 1:n)  {
		names(dims)[i] <- nc$dim[[i]]$name
		if (names(dims)[i]=="sample") {
                  vals <- 1:length(sample.include)
                } else if (names(dims)[i]=="snp") {
                  vals <- snpid[snpsel]
                } else {
                  vals <- nc$dim[[i]]$vals
                }
		dims[[i]] <- dim.def.ncdf(nc$dim[[i]]$name, nc$dim[[i]]$units, vals, unlim=nc$dim[[i]]$unlim)
	}

	# Define variables		
	m <- nc$nvars
	vars <- vector("list",m)
	for(j in 1:m) {
		names(vars)[j] <- nc$var[[j]]$name
		nd <- length(nc$var[[j]]$dim)   	# number of dimensions
		dimsj <- list()
		for(k in 1:nd) {	# get dimensions for jth variable
			dname <- nc$var[[j]]$dim[[k]]$name
			dimsj[[k]] <- dims[[dname]]
		}
		prec <- nc$var[[j]]$prec
		if(prec=="int") prec <- "integer"
		if(prec=="float") prec <- "single"
		vars[[j]] <- var.def.ncdf(nc$var[[j]]$name, nc$var[[j]]$units, dim=dimsj, missval=nc$var[[j]]$missval, prec=prec)
	}
	# Create the netCDF file
	subnc <- create.ncdf(sub.ncdf, vars)
	close.ncdf(subnc)

	# Populate the sub netCDF file
	subnc <- open.ncdf(sub.ncdf, write=TRUE)
		
	# get variable names and their dimensions
	vardim <- list()
	for(j in 1:m){  # for each variable
		vname <- nc$var[[j]]$name
		nd <- length(nc$var[[j]]$dim)   	# number of dimensions
		dname <- rep(NA,nd)
		for(k in 1:nd) dname[k] <- nc$var[[j]]$dim[[k]]$name		# names of dimensions
		vardim[[j]] <- dname
		names(vardim)[j] <- vname
	}

	# for the 1D = "snp" variables, copy data after selecting the snps
	for(j in 1:m) {
		vname <- names(vardim)[j]
		dname <- vardim[[j]]
		chk <- length(dname)==1 & is.element("snp",dname)
		if(!chk) {
                	next
		} else {
			dat <- get.var.ncdf(nc, vname)[snpsel]
			put.var.ncdf(subnc, vname, dat) 
		}
	}

	# for the 1D = "sample" variables, copy data after selecting the samples		
	for(j in 1:m) {
		vname <- names(vardim)[j]
		dname <- vardim[[j]]
		chk <- length(dname)==1 & is.element("sample",dname)
		if(!chk) {
                	next
		} else {
                	dat <- get.var.ncdf(nc, vname)[sampsel]
			put.var.ncdf(subnc, vname, dat) 
		}
	}

	# close to write this to disk and reopen to read "sampleID"
	close.ncdf(subnc)
	subnc <- open.ncdf(sub.ncdf, write=TRUE)
	subid <- get.var.ncdf(subnc, "sampleID")
	ns <- length(subid)

	# for 2D variables, copy one selected sample at a time
	for(j in 1:m) {
		vname <- names(vardim)[j]
		dname <- vardim[[j]]
		chk <- length(dname)==2 & is.element("sample",dname)
		if(!chk) {
                	next
		} else {
			for(k in 1:ns){
				if (verbose & k%%10==0) message(paste(vname, "sample", k, "of", ns))
				index <- which(is.element(sampid, subid[k]))
				dat <- get.var.ncdf(nc, vname, start=c(1,index), count=c(-1,1))[snpsel]				
				put.var.ncdf(subnc, vname, dat, start=c(1,k), count=c(-1,1))
			}
		}
	}
		
	# close and finish
        close.ncdf(nc)
        close.ncdf(subnc)
        return(invisible(NULL))
}




