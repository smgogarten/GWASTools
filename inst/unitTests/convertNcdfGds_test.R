test_convertNcdfGds <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc.orig <- NcdfGenotypeReader(ncfile)
  
  scanID <- getScanID(nc.orig)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- ScanAnnotationDataFrame(data.frame(scanID, sex, stringsAsFactors=FALSE))
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  rsID <- paste("rs", snpID, sep="")
  alleleA <- sample(c("A","C","G","T"), nsnp(nc.orig), replace=TRUE)
  alleleB <- sample(c("A","C","G","T"), nsnp(nc.orig), replace=TRUE)
  snpdf <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position, rsID, alleleA,
                                             alleleB, stringsAsFactors=FALSE))
  geno <- getGenotype(nc.orig)
  # gds stores missing as 3
  geno[is.na(geno)] <- 3
  close(nc.orig)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile, sample.annot=scandf,
                 snp.annot=snpdf)

  gdsobj <- openfn.gds(gdsfile)
  checkEquals(snpID, read.gdsn(index.gdsn(gdsobj, "snp.id")))
  checkEquals(scanID, read.gdsn(index.gdsn(gdsobj, "sample.id")))
  checkEquals(chromosome, read.gdsn(index.gdsn(gdsobj, "snp.chromosome")))
  checkEquals(position, read.gdsn(index.gdsn(gdsobj, "snp.position")))
  checkEquals(rsID, read.gdsn(index.gdsn(gdsobj, "snp.rs.id")))
  checkEquals(paste(alleleA, alleleB, sep="/"),
              read.gdsn(index.gdsn(gdsobj, "snp.allele")))
  checkEquals(geno, read.gdsn(index.gdsn(gdsobj, "genotype")))

  closefn.gds(gdsobj)
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

test_checkNcdfGds <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc.orig <- NcdfGenotypeReader(ncfile)
  
  scanID <- getScanID(nc.orig)
  sex <- c(rep("M", 10), rep("F", 10))
  scandf <- ScanAnnotationDataFrame(data.frame(scanID, sex, stringsAsFactors=FALSE))
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  snpdf <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position))
  geno <- getGenotype(nc.orig)
  # gds stores missing as 3
  geno[is.na(geno)] <- 3
  close(nc.orig)

  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile, sample.annot=scandf,
                 snp.annot=snpdf)

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
