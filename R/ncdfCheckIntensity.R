ncdfCheckIntensity <- function(path=".",
                             intenpath=".",
                             ncdf.filename, 
                             snp.annotation, 
                             scan.annotation, 
                             sep.type,
                             skip.num,
                             col.total,
                             col.nums,
                             scan.name.in.file, 
                             check.scan.index,
                             n.scans.loaded,
                             affy.inten = FALSE,
                             diagnostics.filename = "ncdfCheckIntensity.diagnostics.RData",
                             verbose = TRUE) {

# v4 - add ballelefreq and logrratio and argument N, fix rm snp

	# sx is vector of sample indices to check
	# N is the number of samples loaded so far

	if(!all(is.element(check.scan.index,1:n.scans.loaded))) stop("check.scan.index must be included in 1:n.scans.loaded")
	
	quantfile <- NcdfIntensityReader(ncdf.filename)

	# get sample ids
	quant.sampid <- getVariable(quantfile, "sampleID", start=1, count=n.scans.loaded)

	# get snp ids from the netCDF file
	nc.snpid <- getSnpID(quantfile)
	
	# get sample info and file names
        stopifnot(all(c("scanID", "scanName", "file") %in% names(scan.annotation)))
	
	if(any(!is.element(quant.sampid, scan.annotation$scanID))) stop("some sample id(s) in ncdf file not found in sample annotation dataframe")
	scan.annotation <- scan.annotation[match(quant.sampid, scan.annotation$scanID),]

	if("inten.file" %in% names(scan.annotation)) {
		type <- "Affy"
	} else {
		type <- "Illumina"
	}

	files <- paste(path, scan.annotation$file, sep="/")
	if(type=="Affy") afiles <- paste(intenpath, scan.annotation$inten.file, sep="/")

	# check col.nums vector
        col.nums <- col.nums[!is.na(col.nums)]
        if(!all(names(col.nums) %in% c("snp", "sample", "geno", "a1", "a2", "qs", "x", "y", "rawx", "rawy", "r", "theta", "ballelefreq","logrratio"))) stop("problem with col.nums vector names")
        if(!is.integer(col.nums)) stop("col.nums vector class is not integer")
        if(!("snp" %in% names(col.nums))) stop("snp id missing in col.nums")
        if( max(col.nums) > col.total) stop("some element of col.nums is greater than total number of columns")

	# compare with ncdf
	varin <- names(col.nums)
	vncdf <- getVariableNames(quantfile)
	xn <- c("qs", "x", "y","rawx","rawy","r","theta","ballelefreq","logrratio")
	yn <- c("quality", "X", "Y", "rawX", "rawY", "R", "Theta","BAlleleFreq","LogRRatio")
	zn <- yn[is.element(xn, varin)]
	if(!all(is.element(zn, vncdf))) stop("variables designated in col.nums are not all defined in the ncdf file")
	if(any(is.element(varin, c("geno", "a1", "a2"))) & !is.element("genotype",vncdf)) stop("variables designated in col.nums are not all defined in the ncdf file")

	# check snp.annotation
        stopifnot(all(c("snpID", "snpName") %in% names(snp.annotation)))
	if(any(snp.annotation$snpID != sort(snp.annotation$snpID))) stop("snp annotation ids not in order")
	if(any(snp.annotation$snpID != nc.snpid)) stop("snp annotation ids not the same as in ncdf")
        n <- nrow(snp.annotation)
               
	#generate colClasses vector for read.table
	cc <- rep("NULL",col.total)
	cc[col.nums[names(col.nums) %in% c("snp","sample","geno","a1","a2")]] <- "character"
	cc[col.nums[names(col.nums) %in% c("qs","x","y","rawx","rawy","r","theta","ballelefreq","logrratio")]] <- "double"
	
	#generate names for the genotype data.frame
        df.names <- names(sort(col.nums))

	# set up objects to keep track of things for each file

	# refresh diagnostics from when the ncdf was created
	fn <- length(files)
	read.file <- rep(NA, fn)  # keeps track of whether the file was readable or not
	row.num <- rep(NA, fn)	  # number of rows read
	sample.names <- vector("list",fn)		 # list of vectors of unique sample names in each file
	sample.match <- rep(NA, fn)		# indicator whether sample name inside file matches sample names in sample annotation data.frame
	#missg <- vector("list",fn)		 # vector of character string(s) used for missing genotypes (i.e. not AA, AB or BB)
	snp.chk <- rep(NA,fn)
	chk <- rep(NA,fn)			# final check on data ready to load into ncdf

	# new diagnostics
	snp.order <- rep(NA,fn)
	#geno.chk <- rep(NA,fn)
	qs.chk <- rep(NA,fn)
	z <- rep(NA,fn)
	inten.chk <- list(z,z,z,z,z,z,z,z)
	names(inten.chk) <- c("x", "y", "rawx", "rawy", "r", "theta","ballelefreq","logrratio")	

	# diagnostics for Affy intensity files
	read.file.inten <- rep(NA,fn)
	sample.match.inten <- rep(NA,fn)
	rows.equal <- rep(NA,fn)
	snp.chk.inten <- rep(NA,fn)

	# Set tolerance for intensity comparisons
	tol <- 1e-4 # - used for difference
	rtol <- 1e-6  # used for ratio of difference to raw data file value

	# specify intensity variables to check
	qvars <- c("X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio")
	qvard <- c("x", "y", "rawx", "rawy", "r","theta", "ballelefreq", "logrratio")
	varin <- names(col.nums)
	if(affy.inten==TRUE) varin <- c(varin, "x","y")
	qvars <- qvars[is.element(qvard, varin)]
	qvard <- qvard[is.element(qvard, varin)]
	m <- length(qvars)

	inten.chk <- inten.chk[is.element(names(inten.chk), qvard)]	

	if(is.null(end)) end <- nscan(quantfile)

	nsx <- length(check.scan.index)

        if (verbose) start <- Sys.time()
	for(i in check.scan.index){


		# save diagnostics at each step in case of crash
		diagnostics <- list(read.file, row.num, sample.names, sample.match, snp.chk, chk, snp.order, qs.chk,
        		read.file.inten, sample.match.inten, rows.equal, snp.chk.inten, inten.chk)
		names(diagnostics) <- c("read.file", "row.num", "sample.names", "sample.match", "snp.chk", "chk",
					 "snp.order", "qs.chk", 
					"read.file.inten", "sample.match.inten", "rows.equal", "snp.chk.inten", "inten.chk")
		save(diagnostics, file=diagnostics.filename)

		#message(paste("file",i))

		#read in the file for one sample and keep columns of interest; skip to next file if there is a read error (using function "try")
		if(scan.name.in.file==-1) {skip.num <- skip.num-1; head<-TRUE} else  {head<-FALSE}
		dat <- try(read.table(files[i], header=head, sep=sep.type, comment.char="", skip=skip.num, colClasses=cc))
		if (inherits(dat, "try-error")) { read.file[i] <- 0; message(paste("error reading file",i)); next; k <- k+1 } 		
		read.file[i] <- 1 
		# get sample name from column heading for Affy
		if(scan.name.in.file==-1) {tmp.names <- names(dat)}
		names(dat) <- df.names

		#check and save row number
		row.num[i] <- dim(dat)[1]
		if(row.num[i]!=n) {rm(dat); next; k <- k+1}  # each file should have the same number of rows (one per snp)

		# Sample names for Illumina			
		if(is.element("sample", names(dat))){
			sample.names[[i]] <- unique(dat$sample)
			if(length(sample.names[[i]])>1) {rm(dat);next; k <- k+1}	# there should only be one sample per file
			if(sample.names[[i]]!=scan.annotation$scanName[i]) {sample.match[i] <- 0; rm(dat); next; k <- k+1}  else {sample.match[i] <- 1}
				# sample name inside file should match sample.name vector
		}
		# Sample names for Affy
		if(scan.name.in.file==-1) {
			tmp <- paste(scan.annotation$scanName[i], c("_Call", "_Confidence",".cel"),sep="")
			if(!any(is.element(tmp, tmp.names))) {sample.match[i] <- 0; rm(dat); next} else {sample.match[i] <- 1}
		}	# sample names embedded in file and column names should match

		#check for duplicate snp names
		if(any(duplicated(dat$snp))) {snp.chk[i] <- 0; rm(dat); next; k <- k+1} 

		#check that all expected snps are present
		if(any(!is.element(dat$snp,snp.annotation$snpName))) {snp.chk[i] <- 0; rm(dat); next; k <- k+1} else snp.chk[i] <- 1

		#Using the first raw data file to make it this far, put the int.ids in same order as in raw data
		#	(expecting all to be in this order)
		if(!exists("snp2")) {
			row.names(snp.annotation) <- snp.annotation$snpName
			snp2 <- snp.annotation[dat$snp, ]
		}

		# check to be sure snp ids are in the same order in each file
		if(!all(snp2$snpName==dat$snp)) { rm(dat); snp.order[i] <- 0; next; k <- k+1} else {snp.order[i] <- 1}

		# load quality score from ncdf and check
		if(is.element("qs", names(dat))) {
			qs <- getQuality(quantfile, start=c(1,i), count=c(n,1))
			qs <- qs[match(snp2$snpID, nc.snpid)]
			if(!all(is.na(qs)==is.na(dat$qs)))	{rm(dat); rm(qs); qs.chk[i] <- 0
			} else {
				qs <- qs[!is.na(qs)]
				dqs <- dat$qs[!is.na(dat$qs)]
				dif <- abs(qs-dqs)
				ratio <- dif/dqs
				chkr <- ratio<rtol
				chkd <- dif < tol
				if(all(chkr | chkd)) {qs.chk[i] <- 1; rm(qs); rm(dqs)}  else {rm(dat); rm(qs); rm(dqs); qs.chk[i] <- 0}
			}
		}

		# For Affy, need to replace dat with the intensity data file here 
		if(type=="Affy" & affy.inten==TRUE) {
			rm(dat)

			# Get the intensity data
			dat <- try(read.table(afiles[i], sep="\t", header=TRUE, colClasses=c("character","double")))
			if (inherits(dat, "try-error")) { read.file.inten[i] <- 0; message(paste("error reading intensity file",i)); next }
			read.file.inten[i] <- 1	
	
			# check sample names
				tmp.names <- names(dat)
				tmp <- paste(scan.annotation$scanName[i], c("_Call", "_Confidence",".cel"),sep="")
				if(!any(is.element(tmp, tmp.names))) {sample.match.inten[i] <- 0; rm(dat); next} else {sample.match.inten[i] <- 1}
			names(dat) <- c("sub.probe.id", "inten")
			
			# rearrange to get one row per snp with X and Y for each row
			dat$nchar <- nchar(dat$sub.probe.id)
			dat$sub <- substr(dat$sub.probe.id, dat$nchar, dat$nchar)
			dat$probe.id <- substr(dat$sub.probe.id, 1, dat$nchar-2)
			dat.a <- dat[is.element(dat$sub,"A"),][,c("inten", "probe.id")]; names(dat.a) <- c("x", "probe.id")
			dat.b <- dat[is.element(dat$sub,"B"),][,c("inten", "probe.id")]; names(dat.b) <- c("y", "probe.id")
			if(nrow(dat.a)!=nrow(dat.b) && any(dat.a$probe.id!=dat.b$probe.id)) { rows.equal[i] <- 0; rm(dat); next }
			rows.equal[i] <- 1
			dat <- cbind(dat.a,dat.b[,"y"]); names(dat)[3] <- "y"
			rm(list=(c("dat.a","dat.b")))

			# remove "AFFX snps, check for duplicate snp names, check that all expected snps are present
			dat <- dat[is.element(dat$probe.id,snp2$snpName),]
			if(nrow(dat)!=n) {snp.chk.inten[i] <- 0; rm(dat); next }
			if(any(duplicated(dat$probe.id))) {snp.chk.inten[i] <- 0; rm(dat); next}
			if(any(!is.element(snp2$snpName,dat$probe.id))) {snp.chk.inten[i] <- 0; rm(dat); next} else snp.chk.inten[i] <- 1
		}				
		
		#message(paste("file", i, "point A"))
		
		# load other intensity variable(s) from ncdf and check against what's in the data file
		if(m>0) {
			for(j in 1:m) {
				z <- getVariable(quantfile, varname=qvars[j], start=c(1,i), count=c(n,1))
				z <- z[match(snp2$snpID, nc.snpid)]  # put in the same order as dat
				dz <- dat[,qvard[j]]
				if(!all(is.na(z)==is.na(dz)))	{rm(dat); rm(z); inten.chk[[j]][i] <- 0
				} else {
					z <- z[!is.na(z)]
					dz <- dz[!is.na(dz)]
					dif <- abs(z-dz)
					ratio <- dif/dz
					chkr <- ratio<rtol
					chkd <- dif < tol
					if(all(chkr | chkd)) {inten.chk[[j]][i] <- 1; rm(z); rm(dz)}  else {rm(dat); rm(z); rm(dz); inten.chk[[j]][i] <- 0}
				}
			#message(paste("file", i, "point B", j))
			}
		}

		#message(paste("file", i, "point C"))
		chk[i] <- 1	# made it this far
		rm(dat)
		# to monitor progress
		if(verbose & i%%10==0) {
			rate <- (Sys.time()-start)/10
			percent <- 100*i/nsx
			message(paste("file", i, "-", format(percent,digits=3), "percent completed - rate =", format(rate,digits=4)))
                        start <- Sys.time()
		}

	}	# end of loop




	close(quantfile)
	diagnostics <- list(read.file, row.num, sample.names, sample.match, snp.chk, chk, snp.order, qs.chk,
        		read.file.inten, sample.match.inten, rows.equal, snp.chk.inten, inten.chk)
	names(diagnostics) <- c("read.file", "row.num", "sample.names", "sample.match", "snp.chk", "chk",
					 "snp.order", "qs.chk", 
					"read.file.inten", "sample.match.inten", "rows.equal", "snp.chk.inten", "inten.chk")
	save(diagnostics, file=diagnostics.filename)

	return(diagnostics)

}	# end of check.ncdf
			

