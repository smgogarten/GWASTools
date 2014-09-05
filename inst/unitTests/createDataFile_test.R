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

test_exceptions <- function() {
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- c(101:109, NA)
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  file <- tempfile()

  checkException(createDataFile(annot, file, variables="foo")) # wrong var
  annot <- data.frame(snp=snpID, c=chrom, p=pos) # wrong names
  checkException(createDataFile(annot, file))
  snpID <- paste("rs", 1:10, sep="") # snpID not an integer
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      stringsAsFactors=FALSE)
  checkException(createDataFile(annot, file))
  snpID <- rep(1L, 10) # snpID not unique
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(createDataFile(annot, file))
  snpID <- 10:1 # snpID not sorted
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(createDataFile(annot, file))

  unlink(file)
}

test_createNcdf <- function() {
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- c(101:109, NA)
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  file <- tempfile()
  vars <- c("X","Y")
  nsamp <- 5L
  arrname <- "Array"
  build <- "Build"

  GWASTools:::.createNcdf(annot, file, variables=vars, n.samples=nsamp,
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
  checkIdentical(c(length(snpID), nsamp), dim(getX(nc)))
  checkIdentical(arrname, getAttribute(nc, "array_name"))
  checkIdentical(build, getAttribute(nc, "genome_build"))
  close(nc)

  unlink(file)
}

test_createGds <- function() {
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- c(101:109, NA)
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  file <- tempfile()
  vars <- "genotype"

  GWASTools:::.createGds(annot, file, variables=vars,
             precision="single", compress="ZIP.max")

  gds <- GdsGenotypeReader(file)
  checkEquals(0, nscan(gds))
  checkEquals(length(snpID), nsnp(gds))
  checkIdentical(c("sample.id", "snp.id", "snp.chromosome", "snp.position", vars),
                 getVariableNames(gds))
  checkIdentical(snpID, getSnpID(gds))
  checkIdentical(chrom, getChromosome(gds))
  checkIdentical(pos, getPosition(gds))
  close(gds)

  unlink(file)
}

test_createDataFile <- function() {
  snpID <- 1:5
  chrom <- rep(1L,5)
  pos <- 101:105
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  file <- tempfile()

  createDataFile(annot, file, file.type="gds")
  gds <- GdsGenotypeReader(file)
  close(gds)

  createDataFile(annot, file, file.type="ncdf")
  nc <- NcdfGenotypeReader(file)
  close(nc)

  unlink(file)
}
