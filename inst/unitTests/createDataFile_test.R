
test_var_exceptions <- function() {
  variables <- "genotype"
  intensity.vars <-  c("quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq","LogRRatio")
  col.total <- 10
  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "foo") # bad name
  checkException(GWASTools:::.checkVars(variables, col.nums, col.total, intensity.vars))

  col.nums <- as.character(c(1,2,12,13)) # bad type
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  checkException(GWASTools:::.checkVars(variables, col.nums, col.total, intensity.vars))

  col.nums <- as.character(c(2,12,13))
  names(col.nums) <- c("sample", "a1", "a2") # snp missing
  checkException(GWASTools:::.checkVars(variables, col.nums, col.total, intensity.vars))

  col.nums <- as.integer(c(1,2,12,25)) # bad column
  names(col.nums) <- c("snp", "sample", "a1", "foo")
  checkException(GWASTools:::.checkVars(variables, col.nums, col.total, intensity.vars))

  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "x", "y") # this is a genotype file
  checkException(GWASTools:::.checkVars(variables, col.nums, col.total, intensity.vars))

  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  col.total <- 3 # bad col total
  checkException(GWASTools:::.checkVars(variables, col.nums, col.total, intensity.vars))
}

test_snp_exceptions <- function() {
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- c(101:109, NA)
  name <- paste0("rs", 1:10)
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      snpName=name, stringsAsFactors=FALSE)

  annot <- data.frame(snp=snpID, c=chrom, p=pos) # wrong names
  checkException(GWASTools:::.checkSnpAnnotation(annot))
  snpID <- paste("rs", 1:10, sep="") # snpID not an integer
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      stringsAsFactors=FALSE)
  checkException(GWASTools:::.checkSnpAnnotation(annot))
  snpID <- rep(1L, 10) # snpID not unique
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(GWASTools:::.checkSnpAnnotation(annot))
  snpID <- 10:1 # snpID not sorted
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(GWASTools:::.checkSnpAnnotation(annot))
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

  x <- GWASTools:::.createNcdf(annot, file, variables=vars, n.samples=nsamp,
              precision="single", array.name=arrname,
              genome.build=build)
  close.ncdf(x)

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
  rsID <- paste0("rs", snpID)
  alleleA <- rep("T", 10)
  alleleB <- rep("G", 10)
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      snpName=rsID, alleleA, alleleB, stringsAsFactors=FALSE)
  file <- tempfile()
  vars <- "genotype"

  x <- GWASTools:::.createGds(annot, file, variables=vars,
             precision="single", compress="ZIP.max")
  closefn.gds(x)

  gds <- GdsGenotypeReader(file)
  checkEquals(0, nscan(gds))
  checkEquals(length(snpID), nsnp(gds))
  checkIdentical(c("sample.id", "snp.id", "snp.chromosome", "snp.position",
                   "snp.rs.id", "snp.allele", vars),
                 getVariableNames(gds))
  checkIdentical(snpID, getSnpID(gds))
  checkIdentical(chrom, getChromosome(gds))
  checkIdentical(pos, getPosition(gds))
  ## checkIdentical(rsID, getVariable(gds, "snp.rs.id"))
  checkIdentical(alleleA, getAlleleA(gds))
  checkIdentical(alleleB, getAlleleB(gds))
  close(gds)

  unlink(file)
}


test_affy_ncdf <- function() {
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- createDataFile(path, ncfile, file.type="ncdf", variables="genotype",
                       snpAnnot, scanAnnot, sep.type="\t",
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
}


test_illumina_ncdf <- function() {
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "rsID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  diagfile <- tempfile()
  res <- createDataFile(path, ncfile, file.type="ncdf", variables="genotype",
                       snpAnnot, scanAnnot, sep.type=",",
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

  file.remove(diagfile)
  file.remove(ncfile)
}


test_affy_gds <- function() {
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- createDataFile(path, ncfile, file.type="gds", variables="genotype",
                        snpAnnot, scanAnnot, sep.type="\t",
                       skip.num=1, col.total=6, col.nums=col.nums,
                       scan.name.in.file=-1, diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))

  # check
  nc <- GdsGenotypeReader(ncfile)
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
}


test_illumina_gds <- function() {
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "rsID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  diagfile <- tempfile()
  res <- createDataFile(path, ncfile, file.type="gds",  variables="genotype",
                        snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)

  # check
  nc <- GdsGenotypeReader(ncfile)
  origfile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
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
}


test_intensity_ncdf <- function() {
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "rsID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,5,16,17))
  names(col.nums) <- c("snp", "sample", "quality", "X", "Y")
  diagfile <- tempfile()
  res <- createDataFile(path, ncfile, file.type="ncdf", variables=c("quality", "X", "Y"),
                        snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)

  # check
  nc <- NcdfIntensityReader(ncfile)
  origfile <- system.file("extdata", "illumina_qxy.nc", package="GWASdata")
  nc2 <- NcdfIntensityReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkIdentical(getQuality(nc), getQuality(nc2, snp=c(1,-1), scan=c(1,3)))
  checkIdentical(getX(nc), getX(nc2, snp=c(1,-1), scan=c(1,3)))
  checkIdentical(getY(nc), getY(nc2, snp=c(1,-1), scan=c(1,3)))
  close(nc)
  close(nc2)

  file.remove(diagfile)
  file.remove(ncfile)
}


test_intensity_gds <- function() {
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "rsID", "chromosome", "position")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,5,16,17))
  names(col.nums) <- c("snp", "sample", "quality", "X", "Y")
  diagfile <- tempfile()
  res <- createDataFile(path, ncfile, file.type="gds", variables=c("quality", "X", "Y"),
                        snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)

  # check
  nc <- GdsIntensityReader(ncfile)
  origfile <- system.file("extdata", "illumina_qxy.nc", package="GWASdata")
  nc2 <- NcdfIntensityReader(origfile)
  checkIdentical(getSnpID(nc), getSnpID(nc2))
  checkIdentical(getChromosome(nc), getChromosome(nc2))
  checkIdentical(getPosition(nc), getPosition(nc2))
  checkIdentical(getScanID(nc), getScanID(nc2, 1:3))
  checkIdentical(getQuality(nc), getQuality(nc2, snp=c(1,-1), scan=c(1,3)))
  checkEquals(getX(nc), getX(nc2, snp=c(1,-1), scan=c(1,3)))
  checkEquals(getY(nc), getY(nc2, snp=c(1,-1), scan=c(1,3)))
  close(nc)
  close(nc2)

  file.remove(diagfile)
  file.remove(ncfile)
}

