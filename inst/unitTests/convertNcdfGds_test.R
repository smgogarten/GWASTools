test_convertNcdfGds <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=ncfile, file.type="ncdf")
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
                         n.samples=20, filename=ncfile, file.type="ncdf")
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
  gdsfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=gdsfile, file.type="gds")

  ncfile <- tempfile()
  convertGdsNcdf(gdsfile, ncfile)

  nc.orig <- GdsGenotypeReader(gdsfile)
  nc.new <- NcdfGenotypeReader(ncfile)
  checkEquals(getScanID(nc.orig), getScanID(nc.new))
  checkEquals(getSnpID(nc.orig), getSnpID(nc.new))
  checkEquals(getChromosome(nc.orig), getChromosome(nc.new))
  checkEquals(getPosition(nc.orig), getPosition(nc.new))
  checkEquals(getGenotype(nc.orig), getGenotype(nc.new))

  close(nc.orig)
  close(nc.new)
  file.remove(ncfile, gdsfile)
}

test_convertGdsNcdf_transpose <- function() {
  gdsfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=gdsfile, file.type="gds")

  SNPRelate::snpgdsTranspose(gdsfile)
  ncfile <- tempfile()
  convertGdsNcdf(gdsfile, ncfile)

  gds.orig <- GdsGenotypeReader(gdsfile)
  nc.new <- NcdfGenotypeReader(ncfile)
  checkEquals(getGenotype(gds.orig), getGenotype(nc.new))

  close(nc.new)
  close(gds.orig)
  file.remove(ncfile, gdsfile)
}

test_convertGdsNcdf_intensity <- function() {
  gdsfile <- tempfile()
  simulateIntensityMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=gdsfile, file.type="gds")

  ncfile <- tempfile()
  convertGdsNcdf(gdsfile, ncfile)

  nc.orig <- GdsIntensityReader(gdsfile)
  nc.new <- NcdfIntensityReader(ncfile)
  checkEquals(getQuality(nc.orig), getQuality(nc.new))
  checkEquals(getX(nc.orig), getX(nc.new))
  checkEquals(getY(nc.orig), getY(nc.new))

  close(nc.orig)
  close(nc.new)
  file.remove(ncfile, gdsfile)
}

test_checkNcdfGds <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=ncfile, file.type="ncdf")

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
  nc.new <- nc_open(ncfile, write=TRUE)
  geno <- ncvar_get(nc.new, "genotype")
  geno[geno == 1] <- 2
  ncvar_put(nc.new, "genotype", geno)
  nc_close(nc.new)

  checkTrue(!checkNcdfGds(ncfile, gdsfile))

  file.remove(ncfile, gdsfile)
}

test_convertNcdfGds_chromCodes <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=40, n.samples=5,
                         filename=ncfile, file.type="ncdf")
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
