
.open <- function(filename, file.type) {
    if (file.type == "gds") {
        openfn.gds(filename, readonly=FALSE)
    } else if (file.type == "ncdf") {
        open.ncdf(filename, write=TRUE)
    }
}

.varnames <- function(genofile) {
    if (is(genofile, "gds.class")) {
        ls.gdsn(genofile)
    } else if (is(genofile, "ncdf")) {
        names(genofile$var)
    }
}

.snpid <- function(genofile) {
    if (is(genofile, "gds.class")) {
        read.gdsn(index.gdsn(genofile, "snp.id"))
    } else if (is(genofile, "ncdf")) {
        genofile$dim$snp$vals
    }
}

.addDataGds <- function(genofile, dat, sample.id, vars) {
    append.gdsn(index.gdsn(genofile, "sample.id"), val=sample.id)
    for (v in vars) {
        ## set missing code for genotype
        if (v == "genotype") dat[[v]][is.na(dat[[v]])] <- 3
        append.gdsn(index.gdsn(genofile, v), val=dat[[v]])
    }
}

.addDataNcdf <- function(genofile, dat, sample.id, vars, k, n) {
    put.var.ncdf(genofile, "sampleID", vals=sample.id, start=k, count=1)
    for (v in vars) {
        ## set missing code for genotype
        if (v == "genotype") dat[[v]][is.na(dat[[v]])] <- -1
        put.var.ncdf(genofile, v, vals=dat[[v]], start=c(1,k), count=c(n,1))
    }
}

.addData <- function(genofile, dat, sample.id, vars, k, n) {
    if (is(genofile, "gds.class")) {
        .addDataGds(genofile, dat, sample.id, vars)
    } else if (is(genofile, "ncdf")) {
        .addDataNcdf(genofile, dat, sample.id, vars, k, n)
    }
}

.close <- function(genofile) {
    if (is(genofile, "gds.class")) {
        ## finish compression of nodes
        vars <- .varnames(genofile)
        vars <- vars[!grepl("^snp", vars)] # snp nodes already compressed
        for (v in vars) {
            readmode.gdsn(index.gdsn(genofile, v))
        }

        ## close and cleanup
        filename <- genofile$filename
        closefn.gds(genofile)
        cleanup.gds(filename)
    } else if (is(genofile, "ncdf")) {
        close.ncdf(genofile)
    }
}

addSampleData <- function(path=".",
                          filename,
                          file.type=c("gds", "ncdf"),
                          snp.annotation,
                          scan.annotation,
                          sep.type,
                          skip.num,
                          col.total,
                          col.nums,
                          scan.name.in.file,
                          scan.start.index = 1,
                          diagnostics.filename = "addSampleData.diagnostics.RData",
                          verbose = TRUE) {

    ## get file type
    file.type <- match.arg(file.type)

    genofile <- .open(filename, file.type)

    ## check on variables to read and those in the data file

    ## check col.nums
    col.nums <- col.nums[!is.na(col.nums)]
    intensity.vars <-  c("quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq","LogRRatio")
    if(!all(names(col.nums) %in% c("snp", "sample", "geno", "a1", "a2", intensity.vars))) stop("problem with col.nums vector names")
    if(!is.integer(col.nums)) stop("col.nums vector class is not integer")
    if(!("snp" %in% names(col.nums))) stop("snp id missing in col.nums")
    if( max(col.nums) > col.total) stop("some element of col.nums is greater than total number of columns")

    ## compare with data file
    varnames <- .varnames(genofile)
    chk.vars <- intersect(names(col.nums), intensity.vars)
    if (length(chk.vars) > 0) {
        if(!all(is.element(chk.vars, varnames))) stop("variables designated in col.nums are not all defined in the data file")
    }
    if(any(is.element(names(col.nums), c("geno", "a1", "a2"))) & !is.element("genotype",varnames)) stop("variables designated in col.nums are not all defined in the data file")


    ## Input and check the genotypic data

    ## get sample id information
    stopifnot((c("scanID", "scanName", "file") %in% names(scan.annotation)))
    sample.names <- scan.annotation$scanName
    sample.nums <- scan.annotation$scanID

    ## file names
    data.filenames <- file.path(path, scan.annotation$file)
    rm(scan.annotation)

    ## check snp.annotation
    stopifnot(all(c("snpID", "snpName") %in% names(snp.annotation)))
    if(any(snp.annotation$snpID != sort(snp.annotation$snpID))) stop("snp annotation ids not in order")
    if(any(snp.annotation$snpID != .snpid(genofile))) stop("snp annotation ids not the same as in data file")
    snp.names <- snp.annotation$snpName
    n <- nrow(snp.annotation)
    rm(snp.annotation)

    ## generate colClasses vector for read.table
    cc <- rep("NULL",col.total)
    cc[col.nums[names(col.nums) %in% c("snp","sample","geno","a1","a2")]] <- "character"
    cc[col.nums[names(col.nums) %in% intensity.vars]] <- "double"

    ## generate names for the genotype data.frame
    df.names <- names(sort(col.nums))

    ## set up objects to keep track of things for each file
    fn <- length(data.filenames)
    read.file <- rep(NA, fn)  # indicator of whether the file was readable or not
    row.num <- rep(NA, fn)	  # number of rows read
    samples <- vector("list",fn)	 # list of vectors of unique sample names in each file
    sample.match <- rep(NA, fn)	# indicator whether sample name inside file matches sample.names vector
    missg <- vector("list",fn)	 # vector of character string(s) used for missing genotypes (i.e. not AA, AB or BB)
    snp.chk <- rep(NA,fn)		# indicator for incorrect set of snp ids
    chk <- rep(NA,fn)		# indicator for final check on data ready to load into file
    ## when values for the indicators are assigned, a value of 1 means okay and 0 means not okay
    ## a value of NA means that file processing was aborted before that check was made

    ## add data to file one sample (file) at a time

    k <- scan.start.index

    if (verbose) start <- Sys.time()	# to keep track of the rate of file processing
    for(i in 1:fn){

        ## save diagnostics for each sample in case of interruption
        diagnostics <- list(read.file, row.num, samples, sample.match, missg, snp.chk, chk)
        names(diagnostics) <- c("read.file", "row.num", "samples", "sample.match", "missg", "snp.chk", "chk")
        save(diagnostics, file=diagnostics.filename)

        ## read in the file for one sample and keep columns of interest; skip to next file if there is a read error (using function "try")
        if(scan.name.in.file == -1) {skip.num <- skip.num-1; head<-TRUE} else  {head<-FALSE}
        dat <- try(read.table(data.filenames[i], header=head, sep=sep.type, comment.char="", skip=skip.num, colClasses=cc))
        if (inherits(dat, "try-error")) { read.file[i] <- 0; message(paste("error reading file",i)); next }
        read.file[i] <- 1
        ## get sample name from column heading for Affy
        if(scan.name.in.file == -1) {tmp.names <- names(dat)}
        names(dat) <- df.names
        ##check and save row and sample number info
        row.num[i] <- dim(dat)[1]
        if(row.num[i] != n) {rm(dat); next}  # each file should have the same number of rows (one per snp)
        ## sample names for Illumina (i.e. a sample name column)
        if(is.element("sample", names(dat))){
            samples[[i]] <- unique(dat$sample)
            if(length(samples[[i]])>1) {rm(dat);next}	# there should only be one sample per file
            if(samples[[i]] != sample.names[i]) {sample.match[i] <- 0; rm(dat); next}  else {sample.match[i] <- 1}
            ## sample name inside file should match sample.name vector
        }
        ## sample names for Affy (one column label contains the sample name)
        if(scan.name.in.file == -1) {
            tmp <- paste(sample.names[i], c("_Call", "_Confidence",".cel"),sep="")
            if(!any(is.element(tmp, tmp.names))) {sample.match[i] <- 0; rm(dat); next} else {sample.match[i] <- 1}
        }	## sample name embedded in file name and column name within file should match
        ## check for duplicate snp names
        if(any(duplicated(dat$snp))) {snp.chk[i] <- 0; rm(dat); next}
        ## check that all expected snps are present and sort into int.id order
        if(any(!is.element(dat$snp, snp.names))) {snp.chk[i] <- 0; rm(dat); next} else snp.chk[i] <- 1
        dat <- dat[match(snp.names, dat$snp),]

        if(any(is.element(c("geno","a1","a2"),names(dat)))){
            ## make diploid genotypes if necessary
            if(!is.element("geno", names(dat)) && is.element("a1", names(dat)) && is.element("a2", names(dat))) {
                dat$geno <- paste(dat$a1,dat$a2,sep="")
                new.names <- names(dat)[!is.element(names(dat),c("a1","a2"))]
                dat <- dat[,new.names]
            }
            ## get character string(s) for missing genotypes
            missg[[i]] <- unique(dat$geno[!is.element(dat$geno,c("AA","AB","BB"))])
            ## Make genotypes numeric (number of A alleles in the genotype), with missing = NA
            dat$genotype <- NA
            dat[is.element(dat$geno,"AA"), "genotype"] <- 2
            dat[is.element(dat$geno,c("AB","BA")), "genotype"] <- 1
            dat[is.element(dat$geno,"BB"), "genotype"] <- 0
        }

        ## Load data into file
        vars <- intersect(c("genotype", intensity.vars), names(dat))
        if (i == 1) message("adding variables: ", paste(vars, collapse=", "))
        .addData(genofile, dat, sample.nums[i], vars, k, n)

        k <- k+1	# sample dimension indicator
        chk[i] <- 1	# made it this far
        rm(dat)
        ## to monitor progress
        if(verbose & i%%10==0) {
            rate <- (Sys.time()-start)/10
            percent <- 100*i/fn
            message(paste("file", i, "-", format(percent,digits=3), "percent completed - rate =", format(rate,digits=4)))
            start <- Sys.time()
        }
    }	## end of for loop
    ## finish up
    .close(genofile)
    diagnostics <- list(read.file, row.num, samples, sample.match, missg, snp.chk, chk)
    names(diagnostics) <- c("read.file", "row.num", "samples", "sample.match", "missg", "snp.chk", "chk")
    save(diagnostics, file=diagnostics.filename)
    return(diagnostics)
}

