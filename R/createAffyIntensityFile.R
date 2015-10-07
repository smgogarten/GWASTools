createAffyIntensityFile <- function(path=".",
                                    filename,
                                    file.type=c("gds", "ncdf"),
                                    snp.annotation,
                                    scan.annotation, 
                                    precision = "single",
                                    compress = "ZIP.max",
                                    array.name = NULL,
                                    genome.build = NULL,
                                    diagnostics.filename = "createAffyIntensityFile.diagnostics.RData",
                                    verbose = TRUE) {

    ## scan.start.index is the start index for sample ID and 
    ## n.consecutive.scans is the number of consecutive scans for which to load intensity data (default of -1 means load them all)
    ## filename is the file in which to load the intensity data

    ## get file type
    file.type <- match.arg(file.type)
   
    ## checks
    .checkSnpAnnotation(snp.annotation)
    stopifnot((c("scanID", "scanName", "file") %in% names(scan.annotation)))

    ## create data file
    variables <- c("X", "Y")
    if (file.type == "gds") {
        ## don't need n.samples since we will use append later
        genofile <- .createGds(snp.annotation, filename, variables, precision, compress)
    } else if (file.type == "ncdf") {
        genofile <- .createNcdf(snp.annotation, filename, variables, nrow(scan.annotation),
                                 precision, array.name, genome.build)
    }

    ## get sample id information
    sample.names <- scan.annotation$scanName
    sample.nums <- scan.annotation$scanID
    
    ## file names
    files <- file.path(path, scan.annotation$file)
    fn <- length(files)
    rm(scan.annotation)
    
    ## get snp information
    n <- nrow(snp.annotation)
    snp.names <- snp.annotation$snpName
    rm(snp.annotation)

    ## set up objects to keep track of things for each file
    read.file <- rep(NA, fn)  ## keeps track of whether the file was readable or not
    row.num <- rep(NA, fn)	  ## number of rows read
    rows.equal <- rep(NA,fn)	## rows for A and B equal and ordered the same
    sample.match <- rep(NA,fn)
    snp.chk <- rep(NA,fn)	## all snps are present and no duplicates
    chk <- rep(NA,fn)		  ## final check on data ready to load into ncdf

    ## For each raw data file, rearrange and input to ncdf
    ## this file has two rows per probe, one for A and one for B (sub.probe.id = probe_name-A or probe_name-B)
    ## split the data set and make X=intensity of A and Y=intensity of B

    if (verbose) start <- Sys.time()	## to keep track of the rate of file processing
    for(i in 1:fn){
        
        dat <- try(read.table(files[i], sep="\t", header=TRUE, colClasses=c("character","double")))
        if (inherits(dat, "try-error")) { read.file[i] <- 0; message(paste("error reading file",i)); next }
        read.file[i] <- 1			
        ## check sample names
        tmp.names <- names(dat)
        tmp <- paste0(sample.names[i], c("_Call", "_Confidence", ".cel", ".CEL"))
        if(!any(is.element(tmp, tmp.names))) {sample.match[i] <- 0; rm(dat); next} else {sample.match[i] <- 1}
        names(dat) <- c("sub.probe.id", "inten")
        row.num[i] <- nrow(dat)
	
        ## rearrange to get one row per snp with X and Y for each row
        dat$nchar <- nchar(dat$sub.probe.id)
        dat$sub <- substr(dat$sub.probe.id, dat$nchar, dat$nchar)
        dat$probe.id <- substr(dat$sub.probe.id, 1, dat$nchar-2)
        dat.a <- dat[is.element(dat$sub,"A"),][,c("inten", "probe.id")]; names(dat.a) <- c("X", "probe.id")
        dat.b <- dat[is.element(dat$sub,"B"),][,c("inten", "probe.id")]; names(dat.b) <- c("Y", "probe.id")
        if(nrow(dat.a)!=nrow(dat.b) && any(dat.a$probe.id!=dat.b$probe.id)) { rows.equal[i] <- 0; rm(dat); next }
        rows.equal[i] <- 1
        dat <- cbind(dat.a,dat.b[,"Y"]); names(dat)[3] <- "Y"
        rm(list=(c("dat.a","dat.b")))

        ## remove AFFX snps, check for duplicate snp names, check that all expected snps are present, and sort into int.id order
        dat <- dat[is.element(dat$probe.id,snp.names),]
        if(nrow(dat)!=n) {snp.chk[i] <- 0; rm(dat); next }
        if(any(duplicated(dat$probe.id))) {snp.chk[i] <- 0; rm(dat); next}
        if(any(!is.element(snp.names,dat$probe.id))) {snp.chk[i] <- 0; rm(dat); next} else snp.chk[i] <- 1
        dat <- dat[match(snp.names,dat$probe.id),]

        ## Load data into file
        .addData(genofile, variables, dat, sample.nums[i], i)

        chk[i] <- 1	## made it this far
        rm(dat)
        ## to monitor progress
        if(verbose & i%%10==0) {
            rate <- (Sys.time()-start)/10
            percent <- 100*i/fn
            message(paste("file", i, "-", format(percent,digits=3), "percent completed - rate =", format(rate,digits=4)))
            start <- Sys.time()
        }
    }
    
    .close(genofile, verbose=verbose)
    diagnostics <- list(read.file, row.num, rows.equal, sample.match, snp.chk, chk)
    names(diagnostics) <- c("read.file", "row.num", "rows.equal", "sample.match", "snp.chk", "chk")
    save(diagnostics, file=diagnostics.filename)
    return(diagnostics)
}

