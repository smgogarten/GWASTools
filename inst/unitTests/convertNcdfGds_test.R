test_convertNcdfGds <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc.orig <- NcdfGenotypeReader(ncfile)

  scanID <- getScanID(nc.orig)
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  snpName <- paste("rs", snpID, sep="")
  alleleA <- sample(c("A","C","G","T"), nsnp(nc.orig), replace=TRUE)
  alleleB <- sample(c("A","C","G","T"), nsnp(nc.orig), replace=TRUE)
  snpdf <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position, snpName, alleleA,
                                             alleleB, stringsAsFactors=FALSE))
  geno <- getGenotype(nc.orig)
  close(nc.orig)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile, snp.annot=snpdf)

  gds <- GdsGenotypeReader(gdsfile)
  checkEquals(snpID, getSnpID(gds))
  checkEquals(scanID, getScanID(gds))
  checkEquals(chromosome, getChromosome(gds))
  checkEquals(position, getPosition(gds))
  #checkEquals(snpName, getVariable(gds, "snp.rs.id"))
  checkEquals(alleleA, getAlleleA(gds))
  checkEquals(alleleB, getAlleleB(gds))
  checkEquals(geno, getGenotype(gds))

  close(gds)
  file.remove(ncfile, gdsfile)
}

test_convertNcdfGds_intensity <- function() {
  ncfile <- tempfile()
  simulateIntensityMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc.orig <- NcdfIntensityReader(ncfile)

  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  snpdf <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position))
  close(nc.orig)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile, snp.annot=snpdf)

  nc.orig <- NcdfIntensityReader(ncfile)
  gds <- GdsIntensityReader(gdsfile)
  checkEquals(getQuality(nc.orig), getQuality(gds))
  checkEquals(getX(nc.orig), getX(gds))
  checkEquals(getY(nc.orig), getY(gds))

  close(gds)
  close(nc.orig)
  file.remove(ncfile, gdsfile)
}

test_convertGdsNcdf <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile)
  ncfile2 <- tempfile()
  convertGdsNcdf(gdsfile, ncfile2)

  nc.orig <- NcdfGenotypeReader(ncfile)
  nc.new <- NcdfGenotypeReader(ncfile2)
  checkEquals(getScanID(nc.orig), getScanID(nc.new))
  checkEquals(getSnpID(nc.orig), getSnpID(nc.new))
  checkEquals(getChromosome(nc.orig), getChromosome(nc.new))
  checkEquals(getPosition(nc.orig), getPosition(nc.new))
  checkEquals(getGenotype(nc.orig), getGenotype(nc.new))

  close(nc.orig)
  close(nc.new)
  file.remove(ncfile, gdsfile, ncfile2)
}

test_convertGdsNcdf_intensity <- function() {
  ncfile <- tempfile()
  simulateIntensityMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile)
  ncfile2 <- tempfile()
  convertGdsNcdf(gdsfile, ncfile2)

  nc.orig <- NcdfIntensityReader(ncfile)
  nc.new <- NcdfIntensityReader(ncfile2)
  checkEquals(getQuality(nc.orig), getQuality(nc.new))
  checkEquals(getX(nc.orig), getX(nc.new))
  checkEquals(getY(nc.orig), getY(nc.new))

  close(nc.orig)
  close(nc.new)
  file.remove(ncfile, gdsfile, ncfile2)
}

test_checkNcdfGds <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)

  nc.orig <- NcdfGenotypeReader(ncfile)
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  snpdf <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position))
  close(nc.orig)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile, snp.annot=snpdf)

  checkTrue(checkNcdfGds(ncfile, gdsfile))

  # change ncdf file
  nc.new <- open.ncdf(ncfile, write=TRUE)
  geno <- get.var.ncdf(nc.new, "genotype")
  geno[geno == 1] <- 2
  put.var.ncdf(nc.new, "genotype", geno)
  close.ncdf(nc.new)

  checkTrue(!checkNcdfGds(ncfile, gdsfile))

  file.remove(ncfile, gdsfile)
}

test_convertNcdfGds_chromCodes <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=40, n.samples=5,
                         ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)

  # SNP annotation
  snpdf <- data.frame(snpID=getSnpID(nc),
                      chromosome=getChromosome(nc),
                      position=getPosition(nc))
  snpAnnot <- SnpAnnotationDataFrame(snpdf, autosomeCode=1:38L,
                                     XchromCode=39L, YchromCode=40L,
                                     XYchromCode=41L, MchromCode=42L)
  close(nc)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile, snp.annot=snpAnnot)

  if (require(SNPRelate)) {
    gdsobj <- openfn.gds(gdsfile)
    option <- snpgdsOption(gdsobj)
    checkEquals(1, option$autosome.start)
    checkEquals(38, option$autosome.end)
    checkEquals(39, option$chromosome.code$X)
    checkEquals(40, option$chromosome.code$Y)
    checkEquals(41, option$chromosome.code$XY)
    checkEquals(42, option$chromosome.code$M)
    closefn.gds(gdsobj)
  }
  file.remove(ncfile, gdsfile)
}
