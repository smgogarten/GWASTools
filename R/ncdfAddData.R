ncdfAddData <- function(path="", 
                         ncdf.filename, 
                         snp.annotation,
                         scan.annotation, 
                         sep.type, 
                         skip.num, 
                         col.total, 
                         col.nums, 
                         scan.name.in.file,
                         scan.start.index = 1,
                         diagnostics.filename = "ncdfAddData.diagnostics.RData",
                         verbose = TRUE)
{

	genofile <- open.ncdf(ncdf.filename, write=TRUE)
		
	# check on variables to read and those in the ncdf
        
	# check col.nums
        col.nums <- col.nums[!is.na(col.nums)]
        if(!all(names(col.nums) %in% c("snp", "sample", "geno", "a1", "a2", "qs", "x", "y", "rawx", "rawy", "r", "theta", "ballelefreq","logrratio"))) stop("problem with col.nums vector names")
        if(!is.integer(col.nums)) stop("col.nums vector class is not integer")
	if(!("snp" %in% names(col.nums))) stop("snp id missing in col.nums")
	if( max(col.nums) > col.total) stop("some element of col.nums is greater than total number of columns")
        
	# compare with ncdf
        varin <- names(col.nums)
        vncdf <- names(genofile$var)
        xn <- c("qs", "x", "y","rawx","rawy","r","theta","ballelefreq","logrratio")
        yn <- c("quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq","LogRRatio")
        zn <- yn[is.element(xn, varin)]
        if(!all(is.element(zn, vncdf))) stop("variables designated in col.nums are not all defined in the ncdf file")
        if(any(is.element(varin, c("geno", "a1", "a2"))) & !is.element("genotype",vncdf)) stop("variables designated in col.nums are not all defined in the ncdf file")
			
        
	## Input and check the genotypic data

	# get sample id information
        stopifnot((c("scanID", "scanName", "file") %in% names(scan.annotation)))
	sample.names <- scan.annotation$scanName
	sample.nums <- scan.annotation$scanID
	
	# file names
	data.filenames <- file.path(path, scan.annotation$file)
	rm(scan.annotation)
			
	# check snp.annotation
        stopifnot(all(c("snpID", "snpName") %in% names(snp.annotation)))
	if(any(snp.annotation$snpID != sort(snp.annotation$snpID))) stop("snp annotation ids not in order")
	if(any(snp.annotation$snpID != genofile$dim$snp$vals)) stop("snp annotation ids not the same as in ncdf")
        snp.names <- snp.annotation$snpName
        n <- nrow(snp.annotation)
        rm(snp.annotation)
	
	#generate colClasses vector for read.table
	cc <- rep("NULL",col.total)
	cc[col.nums[names(col.nums) %in% c("snp","sample","geno","a1","a2")]] <- "character"
	cc[col.nums[names(col.nums) %in% c("qs","x","y","rawx","rawy","r","theta","ballelefreq","logrratio")]] <- "double"
	
	#generate names for the genotype data.frame
        df.names <- names(sort(col.nums))
		
	# set up objects to keep track of things for each file
	fn <- length(data.filenames)
	read.file <- rep(NA, fn)  # indicator of whether the file was readable or not
	row.num <- rep(NA, fn)	  # number of rows read
	samples <- vector("list",fn)	 # list of vectors of unique sample names in each file
	sample.match <- rep(NA, fn)	# indicator whether sample name inside file matches sample.names vector
	missg <- vector("list",fn)	 # vector of character string(s) used for missing genotypes (i.e. not AA, AB or BB)
	snp.chk <- rep(NA,fn)		# indicator for incorrect set of snp ids
	chk <- rep(NA,fn)		# indicator for final check on data ready to load into ncdf
	# when values for the indicators are assigned, a value of 1 means okay and 0 means not okay
	# a value of NA means that file processing was aborted before that check was made

	# add data to ncdf one sample (file) at a time

        k <- scan.start.index
        
        if (verbose) start <- Sys.time()	# to keep track of the rate of file processing
	for(i in 1:fn){
          
		# save diagnostics for each sample in case of interruption
		diagnostics <- list(read.file, row.num, samples, sample.match, missg, snp.chk, chk)
		names(diagnostics) <- c("read.file", "row.num", "samples", "sample.match", "missg", "snp.chk", "chk")
		save(diagnostics, file=diagnostics.filename)

		#read in the file for one sample and keep columns of interest; skip to next file if there is a read error (using function "try")
                if(scan.name.in.file == -1) {skip.num <- skip.num-1; head<-TRUE} else  {head<-FALSE}
                dat <- try(read.table(data.filenames[i], header=head, sep=sep.type, comment.char="", skip=skip.num, colClasses=cc))
                if (inherits(dat, "try-error")) { read.file[i] <- 0; message(paste("error reading file",i)); next } 		
                read.file[i] <- 1 
		# get sample name from column heading for Affy
                if(scan.name.in.file == -1) {tmp.names <- names(dat)}
                names(dat) <- df.names
		#check and save row and sample number info
                row.num[i] <- dim(dat)[1]
                if(row.num[i] != n) {rm(dat); next}  # each file should have the same number of rows (one per snp)
		# sample names for Illumina (i.e. a sample name column)
                if(is.element("sample", names(dat))){
                  samples[[i]] <- unique(dat$sample)
                  if(length(samples[[i]])>1) {rm(dat);next}	# there should only be one sample per file
                  if(samples[[i]] != sample.names[i]) {sample.match[i] <- 0; rm(dat); next}  else {sample.match[i] <- 1}
				# sample name inside file should match sample.name vector
                }
		# sample names for Affy (one column label contains the sample name)
                if(scan.name.in.file == -1) {
                         tmp <- paste(sample.names[i], c("_Call", "_Confidence",".cel"),sep="")
                        if(!any(is.element(tmp, tmp.names))) {sample.match[i] <- 0; rm(dat); next} else {sample.match[i] <- 1}
                }	# sample name embedded in file name and column name within file should match		
		#check for duplicate snp names
                if(any(duplicated(dat$snp))) {snp.chk[i] <- 0; rm(dat); next}
		#check that all expected snps are present and sort into int.id order
                if(any(!is.element(dat$snp, snp.names))) {snp.chk[i] <- 0; rm(dat); next} else snp.chk[i] <- 1
                dat <- dat[match(snp.names, dat$snp),]

		if(any(is.element(c("geno","a1","a2"),names(dat)))){
			#make diploid genotypes if necessary
			if(!is.element("geno", names(dat)) && is.element("a1", names(dat)) && is.element("a2", names(dat))) {
                                dat$geno <- paste(dat$a1,dat$a2,sep="")
                                new.names <- names(dat)[!is.element(names(dat),c("a1","a2"))]
                                dat <- dat[,new.names]
                        }
			#get character string(s) for missing genotypes
                        missg[[i]] <- unique(dat$geno[!is.element(dat$geno,c("AA","AB","BB"))])
			# Make genotypes numeric (number of A alleles in the genotype), with missing = -1
                        dat$num <- -1
                        dat[is.element(dat$geno,"AA"), "num"] <- 2
                        dat[is.element(dat$geno,c("AB","BA")), "num"] <- 1
                        dat[is.element(dat$geno,"BB"), "num"] <- 0
		}

		# Load data into ncdf file
                put.var.ncdf(genofile, "sampleID", vals=sample.nums[i], start=k, count=1)
                if(is.element("num",names(dat))) put.var.ncdf(genofile, "genotype", vals=dat$num, start=c(1,k), count=c(n,1))
                if(is.element("qs",names(dat))) put.var.ncdf(genofile, "quality", vals=dat$qs, start=c(1,k), count=c(n,1))
                if(is.element("x",names(dat))) put.var.ncdf(genofile, "X", vals=dat$x, start=c(1,k), count=c(n,1))
                if(is.element("y",names(dat))) put.var.ncdf(genofile, "Y", vals=dat$y, start=c(1,k), count=c(n,1))
                if(is.element("rawx",names(dat))) put.var.ncdf(genofile, "rawX", vals=dat$rawx, start=c(1,k), count=c(n,1))
                if(is.element("rawy",names(dat))) put.var.ncdf(genofile, "rawY", vals=dat$rawy, start=c(1,k), count=c(n,1))
                if(is.element("r",names(dat))) put.var.ncdf(genofile, "R", vals=dat$r, start=c(1,k), count=c(n,1))
                if(is.element("theta",names(dat))) put.var.ncdf(genofile, "Theta", vals=dat$theta, start=c(1,k), count=c(n,1))
                if(is.element("ballelefreq",names(dat))) put.var.ncdf(genofile, "BAlleleFreq", vals=dat$ballelefreq, start=c(1,k), count=c(n,1))
                if(is.element("logrratio",names(dat))) put.var.ncdf(genofile, "LogRRatio", vals=dat$logrratio, start=c(1,k), count=c(n,1))
		
		k <- k+1	# sample dimension indicator
		chk[i] <- 1	# made it this far
		rm(dat)
		# to monitor progress
		if(verbose & i%%10==0) {
			rate <- (Sys.time()-start)/10
			percent <- 100*i/fn
			message(paste("file", i, "-", format(percent,digits=3), "percent completed - rate =", format(rate,digits=4)))
                        start <- Sys.time()
		}
	}	# end of for loop
	# finish up
        close.ncdf(genofile)
        diagnostics <- list(read.file, row.num, samples, sample.match, missg, snp.chk, chk)
        names(diagnostics) <- c("read.file", "row.num", "samples", "sample.match", "missg", "snp.chk", "chk")
        save(diagnostics, file=diagnostics.filename)
        return(diagnostics)
}

