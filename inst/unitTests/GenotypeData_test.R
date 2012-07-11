test_GenotypeData <- function() {
  # simulate data
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)
  scanID <- getScanID(nc)
  sex <- c(rep("M", 10), rep("F", 10))
  plate <- rep(paste("batch", 1:4, sep=""), 5)
  scandf <- data.frame(scanID=scanID, sex=sex, plate=plate)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  snpID <- getSnpID(nc)
  chrom <- getChromosome(nc)
  pos <- getPosition(nc)
  rsID <- paste("rs", snpID, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)

  # creation with only ncdf
  obj <- GenotypeData(nc)
  checkTrue(!hasSnpAnnotation(obj))
  checkTrue(!hasScanAnnotation(obj))
  # creation with only snpAnnot
  obj <- GenotypeData(nc, snpAnnot=snpAnnot)
  checkTrue(hasSnpAnnotation(obj))
  checkTrue(!hasScanAnnotation(obj))
  # creation with only scanAnnot
  obj <- GenotypeData(nc, scanAnnot=scanAnnot)
  checkTrue(!hasSnpAnnotation(obj))
  checkTrue(hasScanAnnotation(obj))
  # creation with both annotations
  obj <- GenotypeData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  checkTrue(hasSnpAnnotation(obj))
  checkTrue(hasScanAnnotation(obj))

  # required variables
  geno <- getGenotype(obj)
  checkIdentical(c(nsnp(obj),nscan(obj)), dim(geno))
  checkIdentical(length(getSnpID(obj)), nsnp(obj))
  checkIdentical(length(getChromosome(obj)), nsnp(obj))
  checkIdentical(length(getPosition(obj)), nsnp(obj))
  checkIdentical(length(getScanID(obj)), nscan(obj))

  # annotation variables - snp
  checkTrue(is.character(getSnpVariableNames(obj)))
  checkTrue(hasSnpVariable(obj, "rsID"))
  checkTrue(!hasSnpVariable(obj, "foo"))
  checkIdentical(NULL, getSnpVariable(obj, "foo"))
  # annotation variables - scan
  checkTrue(is.character(getScanVariableNames(obj)))
  checkTrue(hasSex(obj))
  checkTrue(hasScanVariable(obj, "plate"))
  checkTrue(!hasScanVariable(obj, "foo"))
  checkIdentical(NULL, getScanVariable(obj, "foo"))

  # annotation mismatch
  snpAnnot <- snpAnnot[1:100,]
  checkException(GenotypeData(nc, snpAnnot=snpAnnot))
  scanAnnot$scanID[1] <- 25L
  checkException(GenotypeData(nc, scanAnnot=scanAnnot))

  close(obj)
  unlink(ncfile)
}

test_GenotypeData_Gds <- function() {
  file <- tempfile()
  gds <- createfn.gds(file)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  samp <- 1231:1235
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  add.gdsn(gds, "snp.chromosome", chrom)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "sample.id", samp)
  add.gdsn(gds, "genotype", geno, storage="bit2")
  closefn.gds(gds)

  snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID=snp, chromosome=chrom, position=pos))
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=samp))
  gds <- GdsGenotypeReader(file)
  obj <- GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  geno[geno == 3] <- NA
  checkIdentical(geno, getGenotype(obj))
  checkIdentical(snp, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(samp, getScanID(obj))

  close(obj)
  unlink(file)
}

test_GenotypeData_Matrix <- function() {
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  samp <- 1231:1235
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  mgr <- MatrixGenotypeReader(genotype=geno, snpID=snp, chromosome=chrom, position=pos, scanID=samp)

  snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID=snp, chromosome=chrom, position=pos))
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=samp))
  obj <- GenotypeData(mgr, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  checkIdentical(geno, getGenotype(obj))
  checkIdentical(snp, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(samp, getScanID(obj))
}
