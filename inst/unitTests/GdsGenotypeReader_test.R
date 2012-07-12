test_GdsGenotypeReader <- function() {  
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
  
  obj <- GdsGenotypeReader(file)
  checkIdentical(snp, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(samp, getScanID(obj))
  geno[geno == 3] <- NA
  checkIdentical(geno, getGenotype(obj))

  sel <- samp %in% sample(samp, 3)
  checkIdentical(samp[sel], getScanID(obj, sel))

  chromChar <- getChromosome(obj, char=TRUE)
  checkTrue(is.character(chromChar))
  checkTrue(all(chromChar %in% c(1:22,"X","Y","XY","M","U")))
  checkIdentical(rep(c(1:22,"X","XY","Y","M"), each=10), chromChar)
  close(obj)

  # check alternate chromosome codes
  obj <- GdsGenotypeReader(file, YchromCode=24L, XYchromCode=25L)
  chromChar <- getChromosome(obj, char=TRUE)
  checkIdentical(rep(c(1:22,"X","Y","XY","M"), each=10), chromChar)
  close(obj)

  # check exception with incorrect dimensions
  gds <- openfn.gds(file, readonly=FALSE)
  write.gdsn(index.gdsn(gds, "snp.id"), 1:10)
  closefn.gds(gds)
  checkException(GdsGenotypeReader(file))
  unlink(file)
  
  # check exception with incorrect variable names
  file <- tempfile()
  gds <- createfn.gds(file)
  add.gdsn(gds, "snp.id", snp)
  closefn.gds(gds)
  checkException(GdsGenotypeReader(file))
  unlink(file)
}
