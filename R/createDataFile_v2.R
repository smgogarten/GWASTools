
..checkSnpAnnotation <- function(snp.annotation) {
    ## check that snp annotation has right columns
    stopifnot(all(c("snpID", "chromosome", "position", "snpName") %in% names(snp.annotation)))

    ## make sure all snp annotation columns are integers
    if (!is(snp.annotation$snpID, "integer")) {
        snp.annotation$snpID <- as.integer(snp.annotation$snpID)
        warning(paste("coerced snpID to type integer"))
    }
    if (!is(snp.annotation$chromosome, "integer")) {
        snp.annotation$chromosome <- as.integer(snp.annotation$chromosome)
        warning(paste("coerced chromosome to type integer"))
    }
    if (!is(snp.annotation$position, "integer")) {
        snp.annotation$position <- as.integer(snp.annotation$position)
        warning(paste("coerced position to type integer"))
    }

    ## make sure snpID is unique
    stopifnot(length(snp.annotation$snpID) == length(unique(snp.annotation$snpID)))

    ## make sure snpID is sorted by chromsome and position
    stopifnot(all(snp.annotation$snpID == sort(snp.annotation$snpID)))
    sorted <- order(snp.annotation$chromosome, snp.annotation$position)
    if (!all(snp.annotation$snpID == snp.annotation$snpID[sorted])) {
        stop("snpID is not sorted by chromosome and position")
    }
}

..createNcdf <- function(snp.annotation, filename, variables, n.samples, precision,
                        array.name, genome.build) {

    ## Create the netCDF file and load snp annotation
    ## Define dimensions
    snpdim <- dim.def.ncdf("snp","count", snp.annotation$snpID)
    sampledim <- dim.def.ncdf("sample","count",1:n.samples, unlim=TRUE)

    ## Define variables
    varlist <- list("sampleID"=var.def.ncdf("sampleID", "id", dim=sampledim, missval=0, prec="integer"),
                    "position"=var.def.ncdf("position", "bases", dim=snpdim, missval=-1, prec="integer"),
                    "chromosome"=var.def.ncdf("chromosome", "id", dim=snpdim, missval=-1, prec="integer"))
    if ("genotype" %in% variables) {
        varlist[["genotype"]] <- var.def.ncdf("genotype", "num_A_alleles", dim=list(snpdim,sampledim), missval=-1, prec="byte")
    }
    for (v in setdiff(variables, "genotype")) {
        units <- ifelse(v == "quality", "score", "intensity")
        varlist[[v]] <- var.def.ncdf(v, units, dim=list(snpdim, sampledim), missval=-9999, prec=precision)
    }

    ## Create the netCDF file
    genofile <- create.ncdf(filename, varlist)

    ## Add position data
    put.var.ncdf(genofile, varlist[["position"]], snp.annotation$position)
    put.var.ncdf(genofile, varlist[["chromosome"]], snp.annotation$chromosome)

    ## Add attributes
    if (!is.null(array.name)) att.put.ncdf( genofile, 0, "array_name", array.name )
    if (!is.null(genome.build)) att.put.ncdf( genofile, 0, "genome_build", genome.build )

    sync.ncdf(genofile)
    genofile
}

..createGds <- function(snp.annotation, filename, variables, precision,
                       compress) {

    ## define precision for gds
    precision <- ifelse(precision == "double", "float64", "float32")

    ## define dimensions
    n.snps <- nrow(snp.annotation)

    ## create gds file
    gds <- createfn.gds(filename)

    ## add standard variables
    add.gdsn(gds, "sample.id", storage="integer", valdim=0, compress="") # use valdim=0 then append
    add.gdsn(gds, "snp.id", snp.annotation$snpID, compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.chromosome", snp.annotation$chromosome, storage="uint8",
             compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.position", snp.annotation$position, compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.rs.id", snp.annotation$snpName, compress=compress, closezip=TRUE)
    sync.gds(gds)

    ## add selected variables
    if ("genotype" %in% variables) {
        if (all(c("alleleA", "alleleB") %in% names(snp.annotation))) {
            add.gdsn(gds, "snp.allele", paste(snp.annotation$alleleA, snp.annotation$alleleB, sep="/"),
                     compress=compress, closezip=TRUE)
        }
        geno.node <- add.gdsn(gds, "genotype", storage="bit2", valdim=c(n.snps, 0))
        put.attr.gdsn(geno.node, "snp.order")
    }
    for (v in setdiff(variables, "genotype")) {
        add.gdsn(gds, v, storage=precision, valdim=c(n.snps, 0), compress="")
    }

    sync.gds(gds)
    gds
}

..addData <- function(x, ...) UseMethod(".addData", x)
..addData.gds.class <- function(x, dat, sample.id, vars, ...) {
    append.gdsn(index.gdsn(x, "sample.id"), val=sample.id)
    for (v in vars) {
        ## set missing code for genotype
        if (v == "genotype") dat[[v]][is.na(dat[[v]])] <- 3
        append.gdsn(index.gdsn(x, v), val=dat[[v]])
    }
}

..addData.ncdf <- function(x, dat, sample.id, vars, k, n) {
    put.var.ncdf(x, "sampleID", vals=sample.id, start=k, count=1)
    for (v in vars) {
        ## set missing code for genotype
        if (v == "genotype") dat[[v]][is.na(dat[[v]])] <- -1
        put.var.ncdf(x, v, vals=dat[[v]], start=c(1,k), count=c(n,1))
    }
}

..close <- function(x) UseMethod(".close", x)
..close.gds.class <- function(x) {
    vars <- ls.gdsn(x)
    vars <- vars[!grepl("^snp", vars)] # snp nodes already done
    for (v in vars) readmode.gdsn(index.gdsn(x, v))
    sync.gds(x)

    ## close and cleanup
    filename <- x$filename
    closefn.gds(x)
    cleanup.gds(filename)
}
..close.ncdf <- function(x) close.ncdf(x)

..createDataFile <- function(path=".",
                          filename,
                          file.type=c("gds", "ncdf"),
                            variables = "genotype",
                          snp.annotation,
                          scan.annotation,
                          sep.type,
                          skip.num,
                          col.total,
                          col.nums,
                          scan.name.in.file,
                          scan.start.index = 1,
                           precision = "single",
                           compress = "ZIP.max",
                           array.name = NULL,
                           genome.build = NULL,
                          diagnostics.filename = "createDataFile.diagnostics.RData",
                          verbose = TRUE) {

    ## get file type
    file.type <- match.arg(file.type)

    ## checks

    ## check variables
    intensity.vars <-  c("quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq","LogRRatio")
    stopifnot(all(variables %in% c("genotype", intensity.vars)))

    ## check col.nums
    col.nums <- col.nums[!is.na(col.nums)]
    if(!all(names(col.nums) %in% c("snp", "sample", "geno", "a1", "a2", intensity.vars))) stop("problem with col.nums vector names")
    if(!is.integer(col.nums)) stop("col.nums vector class is not integer")
    if(!("snp" %in% names(col.nums))) stop("snp id missing in col.nums")
    if( max(col.nums) > col.total) stop("some element of col.nums is greater than total number of columns")

    ## compare with variables
    chk.vars <- intersect(names(col.nums), intensity.vars)
    if (length(chk.vars) > 0) {
        if(!all(is.element(chk.vars, variables))) stop("mismatch between col.nums and variables")
    }
    if(any(is.element(names(col.nums), c("geno", "a1", "a2"))) & !is.element("genotype",variables)) stop("mismatch between col.nums and variables")

    ## check snp annotation
    ..checkSnpAnnotation(snp.annotation)

    ## check scan annotation
    stopifnot((c("scanID", "scanName", "file") %in% names(scan.annotation)))

    ## create data file
    if (file.type == "gds") {
        ## don't need n.samples since we will use append later
        genofile <- ..createGds(snp.annotation, filename, variables, precision, compress)
    } else if (file.type == "ncdf") {
        genofile <- ..createNcdf(snp.annotation, filename, variables, nrow(scan.annotation),
                                 precision, array.name, genome.build)
    }


    ## Input and check the genotypic data

    ## get sample id information
    sample.names <- scan.annotation$scanName
    sample.nums <- scan.annotation$scanID

    ## file names
    data.filenames <- file.path(path, scan.annotation$file)
    rm(scan.annotation)

    ## get snp information
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
        ..addData(genofile, dat, sample.nums[i], vars, k, n)

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
    ..close(genofile)
    diagnostics <- list(read.file, row.num, samples, sample.match, missg, snp.chk, chk)
    names(diagnostics) <- c("read.file", "row.num", "samples", "sample.match", "missg", "snp.chk", "chk")
    save(diagnostics, file=diagnostics.filename)
    return(diagnostics)
}

