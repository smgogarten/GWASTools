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

convertGdsNcdf <- function(gds.filename, ncdf.filename,
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

	# define dimensions
	snpdim <- dim.def.ncdf("snp", "count", snp.id)
	sampledim <- dim.def.ncdf("sample", "count", sample.id, unlim=TRUE)

	# define variables
	varID <- var.def.ncdf("sampleID", "id", dim=sampledim, missval=0, prec="integer")
	varpos <- var.def.ncdf("position", "bases", dim=snpdim, missval=-1, prec="integer")
	varchr <- var.def.ncdf("chromosome", "id", dim=snpdim, missval=-1, prec="integer")
	vargeno <- var.def.ncdf("genotype", "num_A_alleles", dim=list(snpdim,sampledim), missval=-1, prec="byte")

	# create the NetCDF file
	ncfile <- create.ncdf(ncdf.filename, list(varID, varpos, varchr, vargeno))

	# add variable data
	if (verbose) message(date(), "\t\twriting position and chromosome ...\n")
	put.var.ncdf(ncfile, varpos, read.gdsn(index.gdsn(gdsobj, "snp.position")))
	put.var.ncdf(ncfile, varchr, read.gdsn(index.gdsn(gdsobj, "snp.chromosome")))

	# add genotype data
	if (verbose) message(date(), "\t\twriting genotypes ...\n")
	node.geno <- index.gdsn(gdsobj, "genotype")
	ndim <- objdesp.gdsn(index.gdsn(gdsobj, "genotype"))$dim
	for (i in 1:ndim[2])
	{
		geno <- read.gdsn(node.geno, start=c(1, i), count=c(-1, 1))
		geno[!(geno %in% c(0,1,2))] <- -1
		put.var.ncdf(ncfile, "genotype", geno, start=c(1,i), count=c(-1,1))
		if (verbose & (i %% 500 == 0))
			message(date(), sprintf("\t\twriting %d samples\n", i))
	}

	# add sample id
	put.var.ncdf(ncfile, varID, sample.id)

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
#   sample.annot  --  the annotation of sample
#   snp.annot  --  the annotation of snp
#   zipflag  --  specify the compression flag except genotype, see "add.gdsn"
#   verbose  --  show progress
#

convertNcdfGds <- function(ncdf.filename, gds.filename,
	sample.annot = NULL, snp.annot = NULL,
        zipflag = "ZIP.max", verbose = TRUE)
{
	# check
	stopifnot(is.character(ncdf.filename))
	stopifnot(is.character(gds.filename))
        if (!is.null(snp.annot)) {
          nc <- NcdfGenotypeReader(ncdf.filename,
                                   autosomeCode=autosomeCode(snp.annot),
                                   XchromCode=XchromCode(snp.annot),
                                   XYchromCode=XYchromCode(snp.annot),
                                   YchromCode=YchromCode(snp.annot),
                                   MchromCode=MchromCode(snp.annot))
        } else {
          nc <- NcdfGenotypeReader(ncdf.filename)
        }
        genoData <- GenotypeData(nc, scanAnnot=sample.annot, snpAnnot=snp.annot)
        close(genoData)

	# start
	if (verbose) message(date(), "\tbegin convertNcdfGds ...\n")

	# open netCDF
	nc <- open.ncdf(ncdf.filename)

	# create GDS file
	gfile <- createfn.gds(gds.filename)

	# the order of samples
	sample.order <- 1:nc$dim$sample$len

	# add "sample.id"
	add.gdsn(gfile, "sample.id", get.var.ncdf(nc, "sampleID")[sample.order],
		compress=zipflag, closezip=TRUE)
	# add "snp.id"
	add.gdsn(gfile, "snp.id", nc$dim$snp$vals, compress=zipflag, closezip=TRUE)
	# add "snp.rs.id"
	if (!is.null(snp.annot) & !is.null(snp.annot[["rsID"]]))
	{
		add.gdsn(gfile, "snp.rs.id", snp.annot[["rsID"]],
			compress=zipflag, closezip=TRUE)
	}
        # add "snp.allele"
	if (!is.null(snp.annot))
	{
                if (!is.null(getAlleleA(snp.annot)) & !is.null(getAlleleB(snp.annot)))
                {
                         allele <- paste(getAlleleA(snp.annot),
                                         getAlleleB(snp.annot), sep="/")
                         add.gdsn(gfile, "snp.allele", allele,
                                  compress=zipflag, closezip=TRUE)
                }
	}
	# add "snp.position"
	add.gdsn(gfile, "snp.position", get.var.ncdf(nc, "position"),
		compress=zipflag, closezip=TRUE)
	# add "snp.chromosome"
	add.gdsn(gfile, "snp.chromosome", get.var.ncdf(nc, "chromosome"), storage="uint8",
		compress=zipflag, closezip=TRUE)
        # add chromosome codes
	if (!is.null(snp.annot))
	{
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
		message(date(), "\tstore sampleID, position, and chromosome.\n")

	# add "genotype", 2 bits to store one genotype
	nSNP <- nc$dim$snp$len
	gGeno <- add.gdsn(gfile, "genotype", storage="bit2",
		valdim=c(nSNP, length(sample.order)))
	put.attr.gdsn(gGeno, "snp.order")

	j <- 1
	for (i in sample.order)
	{
		v <- get.var.ncdf(nc, "genotype", start=c(1, i), count=c(-1, 1))
		v[!(v %in% c(0,1,2))] <- 3
		write.gdsn(gGeno, v, start=c(1, j), count=c(-1, 1))
		if (verbose & (j %% 50 == 0))
			message(date(), "\t", j, "\n")
		j <- j + 1
	}

	# sync file
	sync.gds(gfile)

	# close files
	closefn.gds(gfile)
	close.ncdf(nc)
        cleanup.gds(gds.filename)

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

