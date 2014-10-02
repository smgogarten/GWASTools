test_ncdfCreate <- function() {
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- c(101:109, NA)
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  file <- tempfile()
  vars <- c("X","Y")
  nsamp <- 5
  arrname <- "Array"
  build <- "Build"

  ncdfCreate(annot, file, variables=vars, n.samples=nsamp,
                  precision="single", array.name=arrname,
                  genome.build=build)

  nc <- NcdfIntensityReader(file)
  checkEquals(nsamp, nscan(nc))
  checkEquals(length(snpID), nsnp(nc))
  checkIdentical(c("sampleID", "position", "chromosome", vars),
                 getVariableNames(nc))
  checkIdentical(snpID, getSnpID(nc))
  checkIdentical(chrom, getChromosome(nc))
  checkIdentical(pos, getPosition(nc))
  checkIdentical(arrname, getAttribute(nc, "array_name"))
  checkIdentical(build, getAttribute(nc, "genome_build"))
  close(nc)

  # errors
  checkException(ncdfCreate(annot, file, variables="foo")) # wrong var
  annot <- data.frame(snp=snpID, c=chrom, p=pos) # wrong names
  checkException(ncdfCreate(annot, file))
  snpID <- paste("rs", 1:10, sep="") # snpID not an integer
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      stringsAsFactors=FALSE)
  checkException(ncdfCreate(annot, file))
  snpID <- rep(1L, 10) # snpID not unique
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(ncdfCreate(annot, file))
  snpID <- 10:1 # snpID not sorted
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(ncdfCreate(annot, file))

  unlink(file)
}

test_ncdfAddData <- function() {
  #############
  # Affymetrix
  #############
  # first create empty netCDF
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  ncdfCreate(snpAnnot, ncfile, variables="genotype",
                  n.samples=nrow(scanAnnot))

  # add data
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                       skip.num=1, col.total=6, col.nums=col.nums,
                       scan.name.in.file=-1, diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))

  # check
  nc <- NcdfGenotypeReader(ncfile)
  origfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
  nc2 <- NcdfGenotypeReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkIdentical(getGenotype(nc), getGenotype(nc2, snp=c(1,-1), scan=c(1,3)))
  close(nc)
  close(nc2)

  file.remove(diagfile)
  file.remove(ncfile)


  #############
  # Illumina
  #############
  # first create empty netCDF
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  ncdfCreate(snpAnnot, ncfile, variables="genotype",
                  n.samples=nrow(scanAnnot))

  # add data
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "rsID")]
  names(snpAnnot) <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  diagfile <- tempfile()
  res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)

  # check
  nc <- NcdfGenotypeReader(ncfile)
  origfile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
  nc2 <- NcdfGenotypeReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkIdentical(getGenotype(nc), getGenotype(nc2, snp=c(1,-1), scan=c(1,3)))
  close(nc)
  close(nc2)

  # check error conditions
  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "foo") # bad name
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.character(c(1,2,12,13)) # bad type
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.character(c(2,12,13))
  names(col.nums) <- c("sample", "a1", "a2") # snp missing
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.integer(c(1,2,12,25)) # bad column
  names(col.nums) <- c("snp", "sample", "a1", "foo")
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "x", "y") # this is a genotype file
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  names(snpAnnot) <- c("snpID", "foo") # bad snp names
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  names(snpAnnot) <- c("snpID", "snpName")
  names(scanAnnot) <- c("scanID", "foo", "file") # bad scan names
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  names(scanAnnot) <- c("scanID", "scanName", "file")
  snpAnnot$snpID[1] <- snpAnnot$snpID[2] # bad snp IDs
  checkException({
    res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  file.remove(diagfile)
  file.remove(ncfile)

}

test_ncdfCheckGenotype <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")

  # netCDF
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- ncdfCheckGenotype(path, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1:3,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.order == 1))
  checkTrue(all(res$geno.chk == 1))
  
  # now change a file and see if it detects the difference
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,3] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- ncdfCheckGenotype(tmppath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$geno.chk[1] == 0)

  write.table(orig[1:100,], tmp, sep="\t", row.names=FALSE, quote=FALSE)
  res <- ncdfCheckGenotype(tmppath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$row.num[1] == 100)

  file.remove(tmp,diagfile)
}

test_ncdfCheckIntensity <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file", "inten.file")

  # netCDF
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  intenpath <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  col.nums <- as.integer(c(2,4)); names(col.nums) <- c("snp", "qs")
  diagfile <- tempfile()
  res <- ncdfCheckIntensity(path, intenpath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1:3,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.order == 1))
  checkTrue(all(res$qs.chk == 1))
  for (element in res$inten.chk) checkTrue(all(element == 1))

  # now change a file and see if it detects the difference
  # quality
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,4] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- ncdfCheckIntensity(tmppath, intenpath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  checkTrue(res$qs.chk[1] == 0)

  # intensity
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file", "inten.file")
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$inten.file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,2] <- 0
  tmpfile2 <- "tmp2.txt"
  tmp2 <- paste(tmppath, tmpfile2, sep="/")
  write.table(orig, tmp2, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$inten.file[1] <- tmpfile2
  checkException({
    res <- ncdfCheckIntensity(path, tmppath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  })

  file.remove(tmp, tmp2, diagfile)
}

test_ncdfAddIntensity <- function() {
  # first create empty netCDF
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  ncdfCreate(snpAnnot, ncfile, variables=c("quality","X","Y"),
                  n.samples=nrow(scanAnnot))

  # add sampleID and quality
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")
  scanAnnot1 <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot1) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(2,4)); names(col.nums) <- c("snp", "qs")
  diagfile <- tempfile()
  res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot1, sep.type="\t",
                       skip.num=1, col.total=6, col.nums=col.nums,
                       scan.name.in.file=-1, diagnostics.filename=diagfile)
  file.remove(diagfile)

  # add intensity
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  res <- ncdfAddIntensity(path, ncfile, snpAnnot, scanAnnot)
  checkTrue(all(res$chk == 1))

  # check
  nc <- NcdfIntensityReader(ncfile)
  origfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  nc2 <- NcdfIntensityReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkEquals(getX(nc), getX(nc2, snp=c(1,-1), scan=c(1,3)), tolerance=1e-7)
  checkEquals(getY(nc), getY(nc2, snp=c(1,-1), scan=c(1,3)), tolerance=1e-7)
  close(nc)
  close(nc2)
  
  # check error conditions 
  names(snpAnnot) <- c("snpID", "snpName")
  names(scanAnnot) <- c("scanID", "foo", "file") # bad scan names
  checkException({
    res <- ncdfAddIntensity(path, ncfile, snpAnnot, scanAnnot)
  })
  
  names(scanAnnot) <- c("scanID", "scanName", "file")
  snpAnnot$snpID[1] <- snpAnnot$snpID[2] # bad snp IDs
  checkException({
    res <- ncdfAddIntensity(path, ncfile, snpAnnot, scanAnnot)
  })
  
  file.remove(ncfile)
}
