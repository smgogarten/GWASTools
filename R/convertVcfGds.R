#######################################################################
# Convert a VCF (sequence) file to a GDS file (extract SNP data)
#
# INPUT:
#   vcf.filename -- the file name of VCF format
#   gds.filename -- the output gds file
#   nblock -- the number of lines in buffer
#   compress -- the compression method for sample and snp annotations
#   verbose -- show information
#

convertVcfGds <- function(vcf.filename, gds.filename, nblock=1024, compress="ZIP.max",
	verbose=TRUE)
{
	# check
	stopifnot(is.character(vcf.filename))
	stopifnot(is.character(gds.filename))

	# total number of rows and columns
	Cnt <- count.fields(vcf.filename)
	# check
	if (any(Cnt != Cnt[1])) stop("The file has different numbers of columns.")

	line.cnt <- length(Cnt)
	col.cnt <- max(Cnt)
	if (verbose)
	{
		message("Start convertVcfGds ...\n")
		message("\tOpen ", vcf.filename, "\n")
		message("\tScanning ...\n")
	}


	######################################################################
	# Scan VCF file
	#

	# open the vcf file
	opfile <- file(vcf.filename, open="r")

	# read header
	fmtstr <- substring(readLines(opfile, n=1), 3)
	samp.id <- NULL
	while (length(s <- readLines(opfile, n=1)) > 0)
	{
		if (substr(s, 1, 6) == "#CHROM")
		{
			samp.id <- scan(text=s, what=character(0), quiet=TRUE)[-c(1:9)]
			break
		}
	}
	if (is.null(samp.id))
		stop("Error VCF format: invalid sample id!")

	# scan
	xchr.str <- c(1:22, "X", "Y", "XY", "MT", "M", 23, 24, 25, 26)
	xchr <- as.integer(c(1:22, 23, 25, 24, 26, 26, 23, 24, 25, 26))
	chr <- integer(line.cnt); position <- integer(line.cnt)
	snp.rs <- character(line.cnt)
	snp.allele <- character(line.cnt)
	snp.cnt <- 0

	while (length(s <- readLines(opfile, n=nblock)) > 0)
	{
		for (i in 1:length(s))
		{
			ss <- scan(text=s[i], what=character(0), quiet=TRUE, n=5)
			if (all(ss[c(4,5)] %in% c("A", "G", "C", "T")))
			{
				snp.cnt <- snp.cnt + 1
				chr[snp.cnt] <- xchr[match(ss[1], xchr.str)]
				position[snp.cnt] <- as.integer(ss[2])
				snp.rs[snp.cnt] <- ss[3]
				snp.allele[snp.cnt] <- paste(ss[4], ss[5], sep="/")
			}
		}
	}

	# close the file
	close(opfile)

	# trim
	chr <- chr[1:snp.cnt]; position <- position[1:snp.cnt]
	snp.rs <- snp.rs[1:snp.cnt]; snp.allele <- snp.allele[1:snp.cnt]
	nSamp <- length(samp.id); nSNP <- length(chr)
	geno.str <- c("0|0", "0|1", "1|0", "1|1", "0/0", "0/1", "1/0", "1/1")
	geno.code <- as.integer(c(2, 1, 1, 0, 2, 1, 1, 0))


	######################################################################
	# create GDS file
	#
	gfile <- createfn.gds(gds.filename)

	# add "sample.id"
	add.gdsn(gfile, "sample.id", seq.int(1, length(samp.id)), compress=compress, closezip=TRUE)
	# add "sample.name"
	add.gdsn(gfile, "sample.name", samp.id, compress=compress, closezip=TRUE)
	# add "snp.id"
	add.gdsn(gfile, "snp.id", seq.int(1, length(chr)), compress=compress, closezip=TRUE)
	# add "snp.rs.id"
	add.gdsn(gfile, "snp.rs.id", snp.rs, compress=compress, closezip=TRUE)
	# add "snp.position"
	add.gdsn(gfile, "snp.position", position, compress=compress, closezip=TRUE)
	# add "snp.chromosome"
	add.gdsn(gfile, "snp.chromosome", chr, storage="uint8", compress=compress, closezip=TRUE)
	# add "snp.allele"
	add.gdsn(gfile, "snp.allele", snp.allele, compress=compress, closezip=TRUE)

	# sync file
	sync.gds(gfile)

	if (verbose)
	{
		message(date(), "\tstore sample id, snp id, position, and chromosome.\n")
		message(sprintf("\tstart writing: %d samples, %d SNPs ...\n", nSamp, nSNP))
	}

	# add "gonetype", 2 bits to store one genotype
	gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nSNP, nSamp))
	put.attr.gdsn(gGeno, "snp.order")
	# sync file
	sync.gds(gfile)


	# open the vcf file
	opfile <- file(vcf.filename, open="r")
	# read header
	fmtstr <- substring(readLines(opfile, n=1), 3)
	while (length(s <- readLines(opfile, n=1)) > 0)
	{
		if (substr(s, 1, 6) == "#CHROM")
			break
	}
	# scan
	snp.cnt <- 0
	while (length(s <- readLines(opfile, n=nblock)) > 0)
	{
		for (i in 1:length(s))
		{
			ss <- scan(text=s[i], what=character(0), quiet=TRUE, n=5)
			if (all(ss[c(4,5)] %in% c("A", "G", "C", "T")))
			{
				ss <- scan(text=s[i], what=character(0), quiet=TRUE)[-c(1:9)]
				x <- match(substr(ss, 1, 3), geno.str)
				x <- geno.code[x]
				x[is.na(x)] <- as.integer(3)
				snp.cnt <- snp.cnt + 1
				write.gdsn(gGeno, x, start=c(snp.cnt,1), count=c(1,-1))
			}
		}
	}

	# close the file
	close(opfile)
	# close files
	closefn.gds(gfile)
        cleanup.gds(gds.filename, verbose=verbose)

	if (verbose) message(date(), "\tDone.\n")

	return(invisible(NULL))
}
