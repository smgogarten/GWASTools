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


test_affy_ncdf <- function() {
  # first create empty netCDF
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  createDataFile(snpAnnot, ncfile, file.type="ncdf", variables="genotype",
                  n.samples=nrow(scanAnnot))

  # add data
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
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
  # first create empty netCDF
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  createDataFile(snpAnnot, ncfile, file.type="ncdf", variables="genotype",
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
  res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
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
  # first create empty file
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  createDataFile(snpAnnot, ncfile, file.type="gds", variables="genotype",
                  n.samples=nrow(scanAnnot))

  # add data
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- addSampleData(path, ncfile, file.type="gds", snpAnnot, scanAnnot, sep.type="\t",
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
  # first create empty file
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  createDataFile(snpAnnot, ncfile, file.type="gds", variables="genotype",
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
  res <- addSampleData(path, ncfile, file.type="gds", snpAnnot, scanAnnot, sep.type=",",
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
  # first create empty file
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  createDataFile(snpAnnot, ncfile, file.type="ncdf", variables=c("X", "Y"),
                 n.samples=nrow(scanAnnot))

  # add data
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  snpAnnot <- snpAnnot[,c("snpID", "rsID")]
  names(snpAnnot) <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,16,17))
  names(col.nums) <- c("snp", "sample", "X", "Y")
  diagfile <- tempfile()
  res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
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
  checkIdentical(getX(nc), getX(nc2, snp=c(1,-1), scan=c(1,3)))
  checkIdentical(getY(nc), getY(nc2, snp=c(1,-1), scan=c(1,3)))
  close(nc)
  close(nc2)

  file.remove(diagfile)
  file.remove(ncfile)
}


test_errors <- function() {
  # first create empty netCDF
  data(illumina_snp_annot)
  snpAnnot <- illumina_snp_annot
  data(illumina_scan_annot)
  scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
  ncfile <- tempfile()
  createDataFile(snpAnnot, ncfile, file.type="ncdf", variables="genotype",
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

  # check error conditions
  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "foo") # bad name
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.character(c(1,2,12,13)) # bad type
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.character(c(2,12,13))
  names(col.nums) <- c("sample", "a1", "a2") # snp missing
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.integer(c(1,2,12,25)) # bad column
  names(col.nums) <- c("snp", "sample", "a1", "foo")
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "x", "y") # this is a genotype file
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  col.nums <- as.integer(c(1,2,12,13))
  names(col.nums) <- c("snp", "sample", "a1", "a2")
  names(snpAnnot) <- c("snpID", "foo") # bad snp names
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  names(snpAnnot) <- c("snpID", "snpName")
  names(scanAnnot) <- c("scanID", "foo", "file") # bad scan names
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  names(scanAnnot) <- c("scanID", "scanName", "file")
  snpAnnot$snpID[1] <- snpAnnot$snpID[2] # bad snp IDs
  checkException({
    res <- addSampleData(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type=",",
                       skip.num=11, col.total=21, col.nums=col.nums,
                       scan.name.in.file=1, diagnostics.filename=diagfile)
  })

  file.remove(ncfile)
}
