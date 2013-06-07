ncdfAddIntensity <- function(path=".",
                           ncdf.filename,
                           snp.annotation,
                           scan.annotation, 
                           scan.start.index = 1, 
                           n.consecutive.scans = -1,  
                           diagnostics.filename = "ncdfAddIntensity.diagnostics.RData",
                           verbose = TRUE) 
{

	# scan.start.index is the start index for sample ID and 
	# n.consecutive.scans is the number of consecutive scans for which to load intensity data (default of -1 means load them all)
	# ncdf.filename is the file in which to load the intensity data

		
		
	# open ncdf for writing			
        genofile <- open.ncdf(ncdf.filename, write=TRUE)

	# get sample annotation and names of files to load
	m <- n.consecutive.scans
        sid <- get.var.ncdf(genofile, "sampleID", start=scan.start.index, count=m)
        
        stopifnot(all(c("scanID", "scanName", "file") %in% names(scan.annotation)))
        scan.annotation <- scan.annotation[match(sid, scan.annotation$scanID),]
	sample.names <- scan.annotation$scanName
        files <- file.path(path, scan.annotation$file)
        fn <- length(files)
	rm(scan.annotation)
		
	# check snp.annotation
        stopifnot(all(c("snpID", "snpName") %in% names(snp.annotation)))
	if(any(snp.annotation$snpID != sort(snp.annotation$snpID))) stop("snp annotation ids not in order")
	if(any(snp.annotation$snpID != genofile$dim$snp$vals)) stop("snp annotation ids not the same as in ncdf")
        n <- nrow(snp.annotation)
        snp.names <- snp.annotation$snpName
        rm(snp.annotation)

	# set up objects to keep track of things for each file
        read.file <- rep(NA, fn)  # keeps track of whether the file was readable or not
        row.num <- rep(NA, fn)	  # number of rows read
        rows.equal <- rep(NA,fn)	# rows for A and B equal and ordered the same
        sample.match <- rep(NA,fn)
        snp.chk <- rep(NA,fn)	# all snps are present and no duplicates
        chk <- rep(NA,fn)		  # final check on data ready to load into ncdf

	# For each raw data file, rearrange and input to ncdf
	# this file has two rows per probe, one for A and one for B (sub.probe.id = probe_name-A or probe_name-B)
	# split the data set and make X=intensity of A and Y=intensity of B

        if (verbose) start <- Sys.time()	# to keep track of the rate of file processing
        for(i in 1:fn){
		
                dat <- try(read.table(files[i], sep="\t", header=TRUE, colClasses=c("character","double")))
                if (inherits(dat, "try-error")) { read.file[i] <- 0; message(paste("error reading file",i)); next }
                read.file[i] <- 1			
		# check sample names
                tmp.names <- names(dat)
                tmp <- paste(sample.names[i], c("_Call", "_Confidence",".cel"),sep="")
                if(!any(is.element(tmp, tmp.names))) {sample.match[i] <- 0; rm(dat); next} else {sample.match[i] <- 1}
                names(dat) <- c("sub.probe.id", "inten")
                row.num[i] <- nrow(dat)
	
		# rearrange to get one row per snp with X and Y for each row
                dat$nchar <- nchar(dat$sub.probe.id)
                dat$sub <- substr(dat$sub.probe.id, dat$nchar, dat$nchar)
                dat$probe.id <- substr(dat$sub.probe.id, 1, dat$nchar-2)
                dat.a <- dat[is.element(dat$sub,"A"),][,c("inten", "probe.id")]; names(dat.a) <- c("X", "probe.id")
                dat.b <- dat[is.element(dat$sub,"B"),][,c("inten", "probe.id")]; names(dat.b) <- c("Y", "probe.id")
                if(nrow(dat.a)!=nrow(dat.b) && any(dat.a$probe.id!=dat.b$probe.id)) { rows.equal[i] <- 0; rm(dat); next }
                rows.equal[i] <- 1
                dat <- cbind(dat.a,dat.b[,"Y"]); names(dat)[3] <- "Y"
                rm(list=(c("dat.a","dat.b")))

		# remove AFFX snps, check for duplicate snp names, check that all expected snps are present, and sort into int.id order
                dat <- dat[is.element(dat$probe.id,snp.names),]
                if(nrow(dat)!=n) {snp.chk[i] <- 0; rm(dat); next }
                if(any(duplicated(dat$probe.id))) {snp.chk[i] <- 0; rm(dat); next}
                if(any(!is.element(snp.names,dat$probe.id))) {snp.chk[i] <- 0; rm(dat); next} else snp.chk[i] <- 1
                dat <- dat[match(snp.names,dat$probe.id),]

		# Load data into ncdf file
                j <- scan.start.index + i - 1 		#sample dimension indicator
                put.var.ncdf(genofile, "X", vals=dat$X, start=c(1,j), count=c(n,1))
                put.var.ncdf(genofile, "Y", vals=dat$Y, start=c(1,j), count=c(n,1))

                chk[i] <- 1	# made it this far
                rm(dat)
		# to monitor progress
		if(verbose & i%%10==0) {
			rate <- (Sys.time()-start)/10
			percent <- 100*i/fn
			message(paste("file", i, "-", format(percent,digits=3), "percent completed - rate =", format(rate,digits=4)))
                        start <- Sys.time()
		}
	}
	
        close.ncdf(genofile)
        diagnostics <- list(read.file, row.num, rows.equal, sample.match, snp.chk, chk)
        names(diagnostics) <- c("read.file", "row.num", "rows.equal", "sample.match", "snp.chk", "chk")
        save(diagnostics, file=diagnostics.filename)
        return(diagnostics)
}

