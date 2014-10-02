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

	# add variable data
	if (verbose) message(date(), "\t\tAdding sample data ...\n")
	for (i in 1:length(sample.id)) {
                dat <- list()
                for (v in variables) {
                    node <- index.gdsn(gdsobj, v)
                    dat[[v]] <- read.gdsn(node, start=c(1, i), count=c(-1, 1))
                    if (v == "genotype") dat[[v]][dat[[v]] == 3] <- NA
                }
                .addData(ncfile, dat, sample.id[i], variables, i)
		if (verbose & (i %% 100 == 0))
			message(date(), "\twriting sample ", i, "\n")
	}

	# close files
	close.ncdf(ncfile)
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
#   zipflag  --  specify the compression flag except genotype, see "add.gdsn"
#   verbose  --  show progress
#

convertNcdfGds <- function(ncdf.filename, gds.filename,
	snp.annot = NULL, precision="single",
        zipflag = "ZIP.max", verbose = TRUE)
{
	# check
	stopifnot(is.character(ncdf.filename))
	stopifnot(is.character(gds.filename))

	# start
	if (verbose) message(date(), "\tbegin convertNcdfGds ...\n")

	# open netCDF
	nc <- open.ncdf(ncdf.filename)
        snpID <- nc$dim$snp$vals
        chromosome <- get.var.ncdf(nc, "chromosome")
        position <- get.var.ncdf(nc, "position")
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
        variables <- setdiff(names(nc$var), c("sampleID", "chromosome", "position"))
	gfile <- .createGds(snp.annotation, gds.filename, variables,
                            precision, compress=zipflag)

        # add chromosome codes
	if (!is.null(snp.annot)) {
                put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.start", min(autosomeCode(snp.annot)))
                put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.end", max(autosomeCode(snp.annot)))
                put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "X", XchromCode(snp.annot))
                put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "XY", XYchromCode(snp.annot))
                put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "Y", YchromCode(snp.annot))
                put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "M", MchromCode(snp.annot))
                put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "MT", MchromCode(snp.annot))
	}

	# sync file
	sync.gds(gfile)

	if (verbose)
		message(date(), "\tAdding sample data...\n")

	# add samples
        sample.id <- get.var.ncdf(nc, "sampleID")

	for (i in 1:length(sample.id)) {
                dat <- list()
                for (v in variables) {
                    dat[[v]] <- get.var.ncdf(nc, v, start=c(1, i), count=c(-1, 1))
                }
                .addData(gfile, dat, sample.id[i], variables)
		if (verbose & (i %% 100 == 0))
			message(date(), "\twriting sample ", i, "\n")
	}

	# sync file
	sync.gds(gfile)

	# close files
	closefn.gds(gfile)
	close.ncdf(nc)
        cleanup.gds(gds.filename, verbose=verbose)

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

	nc <- open.ncdf(ncdf.filename)
	gfile <- openfn.gds(gds.filename)

	gdsgeno <- index.gdsn(gfile, "genotype")
	if (is.null(gdsgeno)) {
		message("No genotype exists in the CoreArray GDS file.")
                return(FALSE)
              }
	despgeno <- objdesp.gdsn(gdsgeno)
	if (length(despgeno$dim) != 2) {
		message("The dimension of genotype data in the CoreArray GDS file is not correct.")
                return(FALSE)
              }

	# # of samples
	if (nc$dim$sample$len != despgeno$dim[2]) {
		message("The numbers of samples are not equal!")
                return(FALSE)
              }
	else
		if (verbose) message("The numbers of samples are equal\n")

	# # of snps
	if (nc$dim$snp$len != despgeno$dim[1]) {
		message("The numbers of snps are not equal!")
                return(FALSE)
              }
	else
		if (verbose) message("The numbers of snps are equal.\n")

	# checking block by block
	if (verbose) message(date(), "\tChecking in the sample order ...\n")
	for (i in 1:nc$dim$sample$len)
	{
		gN <- get.var.ncdf(nc, "genotype", start=c(1, i), count=c(-1, 1))
		gG <- read.gdsn(gdsgeno, start=c(1, i), count=c(-1, 1))
		gG[gG==3] <- -1
		if (sum(gN != gG) > 0) {
			message(sprintf("The %dth sample error!", i))
                        return(FALSE)
                      }
		if (i %% 500 == 0)
			if (verbose) message(date(), "\t", i, "OK!\n")
	}

	# checking block by block
	if (verbose) message(date(), "\tChecking in the snp order ...\n")
	breaks <- 5000
	bkNum <- nc$dim$snp$len %/% breaks
	bkEnd <- nc$dim$snp$len %% breaks
	sm <- 0
	if (bkNum > 0)
	{
		for (i in 1:bkNum)
		{
			gN <- get.var.ncdf(nc, "genotype", start=c((i-1)*breaks+1, 1), count=c(breaks, -1))
			gG <- read.gdsn(gdsgeno, start=c((i-1)*breaks+1, 1), count=c(breaks, -1))
			gG[gG==3] <- -1
			if (sum(gN != gG) > 0) {
				message(sprintf("The %dth - %dth snps error!", (i-1)*breaks+1, i*breaks))
                                return(FALSE)
                              }
			if (i%%10 == 0)
				if (verbose) message(date(), "\t", i*breaks, "OK!\n")
		}
	}
	if (bkEnd > 0)
	{
		gN <- get.var.ncdf(nc, "genotype", start=c(bkNum*breaks+1, 1), count=c(bkEnd, -1))
		gG <- read.gdsn(gdsgeno, start=c(bkNum*breaks+1, 1), count=c(bkEnd, -1))
		gG[gG==3] <- -1
		if (sum(gN != gG) > 0) {
			message(sprintf("The %dth - %dth snps error!", bkNum*breaks+1, bkNum*breaks+bkEnd))
                        return(FALSE)
                      }
	}

	if (verbose) message("OK!!!\n")
	closefn.gds(gfile)
	close.ncdf(nc)
	return(TRUE)
}

