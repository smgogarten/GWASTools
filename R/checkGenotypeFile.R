checkGenotypeFile <- function(path=".",
                            filename, 
                            file.type=c("gds", "ncdf"),
                            snp.annotation, 
                            scan.annotation, 
                            sep.type,
                            skip.num,
                            col.total,
                            col.nums,
                            scan.name.in.file, 
                            check.scan.index,
                            n.scans.loaded,
                            diagnostics.filename = "checkGenotypeFile.diagnostics.RData",
                            verbose = TRUE) {
		
	# sx is vector of sample indices to check
	# N is the number of samples loaded so far

	if(!all(is.element(check.scan.index,1:n.scans.loaded))) stop("check.scan.index must be included in 1:n.scans.loaded")

        ## get file type
        file.type <- match.arg(file.type)

        if (file.type == "gds") {
            genofile <- GdsGenotypeReader(filename)
        } else if (file.type == "ncdf") {
            genofile <- NcdfGenotypeReader(filename)
        }

	# get sample and snp ids
	geno.sampid <- getScanID(genofile, index=1:n.scans.loaded)
	nc.snpid <- getSnpID(genofile)
	
	# get sample info and file names
        stopifnot(all(c("scanID", "scanName", "file") %in% names(scan.annotation)))
	
	if(any(!is.element(geno.sampid, scan.annotation$scanID))) stop("some sample id(s) in ncdf file not found in sample annotation dataframe")
	scan.annotation <- scan.annotation[match(geno.sampid, scan.annotation$scanID),]
	files <- file.path(path, scan.annotation$file)

	# check col.nums vector
        col.nums <- col.nums[!is.na(col.nums)]
        intensity.vars <-  c("quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq","LogRRatio")
        if(!all(names(col.nums) %in% c("snp", "sample", "geno", "a1", "a2", intensity.vars))) stop("problem with col.nums vector names")
        if(!is.integer(col.nums)) stop("col.nums vector class is not integer")
        if(!("snp" %in% names(col.nums))) stop("snp id missing in col.nums")
        if( max(col.nums) > col.total) stop("some element of col.nums is greater than total number of columns")
	
	# check snp.annotation
        stopifnot(all(c("snpID", "snpName") %in% names(snp.annotation)))
	if(any(snp.annotation$snpID != sort(snp.annotation$snpID))) stop("snp annotation ids not in order")
	if(any(snp.annotation$snpID != nc.snpid)) stop("snp annotation ids not the same as in ncdf")
        n <- nrow(snp.annotation)
                
	#generate colClasses vector for read.table
	cc <- rep("NULL",col.total)
	cc[col.nums[names(col.nums) %in% c("snp","sample","geno","a1","a2")]] <- "character"
        # don't need these to check genotype! 
	#cc[col.nums[names(col.nums) %in% intensity.vars]] <- "double"  
	
	#generate names for the genotype data.frame
        df.names <- names(sort(col.nums))
		
	# set up objects to keep track of things for each file

	# repeat diagnostics from when the ncdf was created
	fn <- length(files)
	read.file <- rep(NA, fn)  # keeps track of whether the file was readable or not
	row.num <- rep(NA, fn)	  # number of rows read
	sample.names <- vector("list",fn)		 # list of vectors of unique sample names in each file
	sample.match <- rep(NA, fn)		# indicator whether sample name inside file matches sample names in sample annotation data.frame
	missg <- vector("list",fn)		 # vector of character string(s) used for missing genotypes (i.e. not AA, AB or BB)
	snp.chk <- rep(NA,fn)
	chk <- rep(NA,fn)			# final check on data ready to load into ncdf

	# new diagnostics
	## snp.order <- rep(NA,fn)
	geno.chk <- rep(NA,fn)

	if(is.null(end)) end <- nscan(genofile)

	nsx <- length(check.scan.index)

        if (verbose) start <- Sys.time()	# to keep track of the rate of file processing
	for(i in check.scan.index){


		# save at each iteration in case of crash
		diagnostics <- list(read.file, row.num, sample.names, sample.match, missg, snp.chk, chk, geno.chk)
		names(diagnostics) <- c("read.file", "row.num", "sample.names", "sample.match", "missg", "snp.chk", "chk",
					 "geno.chk")
		save(diagnostics, file=diagnostics.filename)

		#read in the file for one sample and keep columns of interest; skip to next file if there is a read error (using function "try")
		if(scan.name.in.file==-1) {skip.num <- skip.num-1; head<-TRUE} else  {head<-FALSE}
		dat <- try(read.table(files[i], header=head, sep=sep.type, comment.char="", skip=skip.num, colClasses=cc))
		if (inherits(dat, "try-error")) { read.file[i] <- 0; message(paste("error reading file",i)); next } 		
		read.file[i] <- 1 
		# get sample name from column heading for Affy
		if(scan.name.in.file==-1) {tmp.names <- names(dat)}
		names(dat) <- df.names

		#check and save row number
		row.num[i] <- dim(dat)[1]
		if(row.num[i]!=n) {rm(dat); next}  # each file should have the same number of rows (one per snp)

		# Sample names for Illumina			
		if(is.element("sample", names(dat))){
			sample.names[[i]] <- unique(dat$sample)
			if(length(sample.names[[i]])>1) {rm(dat);next}	# there should only be one sample per file
			if(sample.names[[i]]!=scan.annotation$scanName[i]) {sample.match[i] <- 0; rm(dat); next}  else {sample.match[i] <- 1}
				# sample name inside file should match sample.name vector
		}
		# Sample names for Affy
		if(scan.name.in.file==-1) {
			tmp <- paste(scan.annotation$scanName[i], c("_Call", "_Confidence",".cel"),sep="")
			if(!any(is.element(tmp, tmp.names))) {sample.match[i] <- 0; rm(dat); next} else {sample.match[i] <- 1}
		}	# sample names embedded in file and column names should match

		#check for duplicate snp names
		if(any(duplicated(dat$snp))) {snp.chk[i] <- 0; rm(dat); next} 

		#check that all expected snps are present
		if(!setequal(dat$snp,snp.annotation$snpName)) {snp.chk[i] <- 0; rm(dat); next} else snp.chk[i] <- 1

		#make diploid genotypes if necessary
                if(!is.element("geno", names(dat)) && is.element("a1", names(dat)) && is.element("a2", names(dat))) {
                        dat$geno <- paste(dat$a1,dat$a2,sep="")
                        new.names <- names(dat)[!is.element(names(dat),c("a1","a2"))]
                        dat <- dat[,new.names]}
		#get character string(s) for missing genotypes - this only works when there is only one code for missing genotype
		#missg[[i]] <- unique(dat$geno[!is.element(dat$geno,c("AA","AB","BB"))])
		#if(length(missg[[i]])!=1) { rm(dat); next }
		#make all missing genotypes blank
                missg[[i]] <- ""
                dat[!is.element(dat$geno,c("AA","AB","BB")),"geno"] <- missg[[i]]

		#Using the first raw data file to make it this far, put the int.ids in same order as in raw data
		#	(expecting all to be in this order)
                dat <- dat[match(snp.annotation$snpName, dat$snp),]
		## if(!exists("snp2")) {
		## 	row.names(snp.annotation) <- snp.annotation$snpName
		## 	snp2 <- snp.annotation[dat$snp, ]
		## }

		# check to be sure snp ids are in the same order in each file
		## if(!all(snp2$snpName==dat$snp)) { rm(dat); snp.order[i] <- 0; next} else {snp.order[i] <- 1}

		#load genotypes from ncdf
		geno <- getGenotype(genofile, snp=c(1,n), scan=c(i,1))

		#put in same order as snps in raw data file
		## geno <- geno[match(snp2$snpID, nc.snpid)]
					
		# convert to AB type
		abtype <- rep(NA, n)
		abtype[is.na(geno)] <- missg[[i]]
		abtype[geno==2] <- "AA"
		abtype[geno==1] <- "AB"
		abtype[geno==0] <- "BB"
		if(length(abtype[is.na(abtype)])!=0) {rm(dat); geno.chk[i] <- 0; next }

		# compare genotypes
		if(all(abtype==dat$geno)){geno.chk[i] <- 1; rm(geno); rm(abtype)
			} else {rm(dat); rm(geno); rm(abtype); geno.chk[i] <- 0; next}

		chk[i] <- 1	# made it this far
		if (exists("dat")) rm(dat)
		# to monitor progress
		if(verbose & i%%10==0) {
			rate <- (Sys.time()-start)/10
			percent <- 100*i/nsx
			message(paste("file", i, "-", format(percent,digits=3), "percent completed - rate =", format(rate,digits=4)))
                        start <- Sys.time()
		}


	}	# end of loop

	diagnostics <- list(read.file, row.num, sample.names, sample.match, missg, snp.chk, chk,geno.chk)
	names(diagnostics) <- c("read.file", "row.num", "sample.names", "sample.match", "missg", "snp.chk", "chk",
					"geno.chk")
	save(diagnostics, file=diagnostics.filename)

	close(genofile)
	return(diagnostics)

}	
			

