#######################################################################
# Author: Xiuwen Zheng
# Email: zhengx@u.washington.edu
#



#######################################################################
# To convert a genotype GDS file to netCDF format
#
# INPUT:
#   gdsobj -- a object of gds file
#   ncdf.filename -- the file name of genotype in netCDF format
#   verbose -- show progress
#

convertGdsNcdf <- function(gds.filename, ncdf.filename, precision = "single",
                           verbose = TRUE)
{
	# check
	stopifnot(is.character(gds.filename))
	stopifnot(is.character(ncdf.filename))

	# start
	if (verbose) message(date(), "\tbegin convertGdsNcdf ...\n")

        # open GDS file
        gdsobj <- openfn.gds(gds.filename)

	# read data from a GDS file
	snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	if (!is.numeric(snp.id)) snp.id <- 1:length(snp.id)
	sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	if (!is.numeric(sample.id)) sample.id <- 1:length(sample.id)

        snp.annot <- data.frame(snpID=snp.id,
                                chromosome=read.gdsn(index.gdsn(gdsobj, "snp.chromosome")),
                                position=read.gdsn(index.gdsn(gdsobj, "snp.position")))

	# create the NetCDF file
	if (verbose) message(date(), "\t\tCreating NetCDF file ...\n")
        variables <- ls.gdsn(gdsobj)
        variables <- variables[!grepl("^sample", variables)]
        variables <- variables[!grepl("^snp", variables)]
	ncfile <- .createNcdf(snp.annot, ncdf.filename, variables,
                              n.samples=length(sample.id),
                              precision, array.name=NULL, genome.build=NULL)

        # decide whether gds is (snp,sample) or (sample,snp)
        if ("genotype" %in% variables) {
            dim.attr <- get.attr.gdsn(index.gdsn(gdsobj, "genotype"))
            if ("snp.order" %in% names(dim.attr)) {
                dim1 <- "snp"
            } else if ("sample.order" %in% names(dim.attr)) {
                dim1 <- "sample"
            } else {
                # look at genotype dimensions
                dim.geno <- objdesp.gdsn(index.gdsn(gdsobj, "genotype"))$dim
                if (all(dim.geno == c(length(snp.id), length(sample.id)))) {
                    dim1 <- "snp"
                } else if (all(dim.geno == c(length(sample.id), length(snp.id)))) {
                    dim1 <- "sample"
                } else {
                    stop("Cannot determine dimensions of genotype node")
                }
            }
        } else {
            dim1 <- "snp"
        }
        
	# add variable data
	if (verbose) message(date(), "\t\tAdding sample data ...\n")
	for (i in 1:length(sample.id)) {
                dat <- list()
                for (v in variables) {
                    node <- index.gdsn(gdsobj, v)
                    if (dim1 == "snp") {
                        start <- c(1,i)
                        count <- c(-1,1)
                    } else {
                        start <- c(i,1)
                        count <- c(1,-1)
                    }
                    dat[[v]] <- read.gdsn(node, start=start, count=count)
                    if (v == "genotype") dat[[v]][dat[[v]] == 3] <- NA
                }
                .addData(ncfile, variables, dat, sample.id[i], i)
		if (verbose & (i %% 100 == 0))
			message(date(), "\twriting sample ", i, "\n")
	}

	# close files
	.close(ncfile)
        closefn.gds(gdsobj)

	if (verbose) message(date(), "\tend convertGdsNcdf.\n")
        return(invisible(NULL))
}


#######################################################################
# To convert a genotype netCDF file to CoreArray format
#
# INPUT:
#   ncdf.filename  --  the input file name of genotype in netCDF format
#   gds.filename  --  the output file name of genotype in CoreArray GDS format
#   snp.annot  --  the annotation of snp
#   compress  --  specify the compression flag except genotype, see "add.gdsn"
#   verbose  --  show progress
#

convertNcdfGds <- function(ncdf.filename, gds.filename,
	snp.annot = NULL, precision="single",
        compress = "ZIP_RA", verbose = TRUE)
{
	# check
	stopifnot(is.character(ncdf.filename))
	stopifnot(is.character(gds.filename))

	# start
	if (verbose) message(date(), "\tbegin convertNcdfGds ...\n")

	# open netCDF
	nc <- NcdfReader(ncdf.filename)
        snpID <- getVariable(nc, "snp")
        chromosome <- getVariable(nc, "chromosome")
        position <- getVariable(nc, "position")
	if (!is.null(snp.annot)) {
                stopifnot(allequal(snp.annot$snpID, snpID))
                stopifnot(allequal(snp.annot$chromosome, chromosome))
                stopifnot(allequal(snp.annot$position, position))
                snp.annotation <- pData(snp.annot)
        } else {
                snp.annotation <- data.frame(snpID, chromosome, position)
        }

	# create GDS file
	if (verbose) message(date(), "\tCreating GDS file ...\n")
        variables <- setdiff(getVariableNames(nc), c("sampleID", "chromosome", "position"))
	gfile <- .createGds(snp.annotation, gds.filename, variables,
                            precision, compress=compress)

        # add chromosome codes
	if (!is.null(snp.annot)) {
            .addChromosomeAttributes(gfile, snp.annot)
	}

	# sync file
	sync.gds(gfile)

	if (verbose)
		message(date(), "\tAdding sample data...\n")

	# add samples
        sample.id <- getVariable(nc, "sampleID")

	for (i in 1:length(sample.id)) {
                dat <- list()
                for (v in variables) {
                    dat[[v]] <- getVariable(nc, v, start=c(1, i), count=c(-1, 1))
                }
                .addData(gfile, variables, dat, sample.id[i])
		if (verbose & (i %% 100 == 0))
			message(date(), "\twriting sample ", i, "\n")
	}

	# sync file
	sync.gds(gfile)

	# close files
	close(nc)
        .close(gfile)

	if (verbose)
		message(date(), "\tend convertNcdfGds.\n")

	return(invisible(NULL))
}



#############################################################################
#
# To check the consistence of genotype comparing a netCDF file to a CoreArray file
#
#    ncdf.filename  --  netCDF file name of genotype
#    gds.filename  --  CoreArray gds file name of genotype
#
#############################################################################

checkNcdfGds <- function(ncdf.filename, gds.filename, verbose = TRUE)
{
	stopifnot(is.character(ncdf.filename))
	stopifnot(is.character(gds.filename))

	if (verbose) message(date(), "\tbegin checkNcdfGds ...\n")

	nc <- NcdfGenotypeReader(ncdf.filename)
	gfile <- GdsGenotypeReader(gds.filename)

	if (!hasVariable(gfile, "genotype")) {
		message("No genotype exists in the CoreArray GDS file.")
                return(FALSE)
              }
	if (length(getDimension(gfile, "genotype")) != 2) {
		message("The dimension of genotype data in the CoreArray GDS file is not correct.")
                return(FALSE)
              }

	# # of samples
        nscan <- length(getScanID(nc))
	if (nscan != length(getScanID(gfile))) {
		message("The numbers of samples are not equal!")
                return(FALSE)
              }
	else
		if (verbose) message("The numbers of samples are equal\n")

	# # of snps
        nsnp <- length(getSnpID(nc))
	if (nsnp != length(getSnpID(gfile))) {
		message("The numbers of snps are not equal!")
                return(FALSE)
              }
	else
		if (verbose) message("The numbers of snps are equal.\n")

	# checking block by block
	if (verbose) message(date(), "\tChecking in the sample order ...\n")
	for (i in 1:nscan)
	{
		gN <- getGenotype(nc, snp=c(1,-1), scan=c(i,1))
		gG <- getGenotype(gfile, snp=c(1,-1), scan=c(i,1))
		if (!allequal(gN, gG)) {
			message(sprintf("The %dth sample error!", i))
                        return(FALSE)
                      }
		if (i %% 500 == 0)
			if (verbose) message(date(), "\t", i, "OK!\n")
	}	

	if (verbose) message("OK!!!\n")
	close(gfile)
	close(nc)
	return(TRUE)
}

