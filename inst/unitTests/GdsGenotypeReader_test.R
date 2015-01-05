test_GdsGenotypeReader <- function() {  
  file <- tempfile()
  gds <- createfn.gds(file)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  a <- rep("A", 260)
  b <- rep("G", 260)
  alleles <- paste(a, b, sep="/")
  samp <- 1231:1235
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  add.gdsn(gds, "snp.chromosome", chrom)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  add.gdsn(gds, "genotype", geno, storage="bit2")
  closefn.gds(gds)
  
  obj <- GdsGenotypeReader(file)
  checkIdentical(obj@genotypeDim, "snp,scan")
  checkIdentical(snp, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(a, getAlleleA(obj))
  checkIdentical(b, getAlleleB(obj))
  checkIdentical(samp, getScanID(obj))
  geno[geno == 3] <- NA
  checkIdentical(geno, getGenotype(obj))
  checkIdentical(t(geno), getGenotype(obj, transpose=TRUE))
  # check a subset
  checkIdentical(geno[1:100, 1:3], getGenotype(obj, snp=c(1,100), scan=c(1,3)))
  checkIdentical(t(geno[1:100, 1:3]), getGenotype(obj, snp=c(1,100), scan=c(1,3), transpose=TRUE))
  
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
  
  # check using an existing gds object
  gds <- openfn.gds(file)
  obj <- GdsGenotypeReader(gds)
  checkIdentical(snp, getSnpID(obj))
  closefn.gds(gds)
  
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

test_genotypeDim <- function(){
  
  file <- tempfile()
  gds <- createfn.gds(file)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  a <- rep("A", 260)
  b <- rep("G", 260)
  alleles <- paste(a, b, sep="/")
  samp <- 1231:1235
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  add.gdsn(gds, "snp.chromosome", chrom)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  add.gdsn(gds, "genotype", t(geno), storage="bit2") # transpose of geno here for genotypeDim=scan,snp
  closefn.gds(gds)
  
  obj <- GdsGenotypeReader(file)
  checkIdentical(obj@genotypeDim, "scan,snp")
  checkIdentical(snp, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(a, getAlleleA(obj))
  checkIdentical(b, getAlleleB(obj))
  checkIdentical(samp, getScanID(obj))
  geno[geno == 3] <- NA
  checkIdentical(geno, getGenotype(obj))
  checkIdentical(t(geno), getGenotype(obj, transpose=TRUE))
  # check a subset
  checkIdentical(geno[1:100, 1:3], getGenotype(obj, snp=c(1,100), scan=c(1,3)))
  checkIdentical(t(geno[1:100, 1:3]), getGenotype(obj, snp=c(1,100), scan=c(1,3), transpose=TRUE))
  
  sel <- samp %in% sample(samp, 3)
  checkIdentical(samp[sel], getScanID(obj, sel))
  close(obj)
  
  checkException(GdsGenotypeReader(file, genotypeDim="notavalidstring"))
  checkException(GdsGenotypeReader(file, genotypeDim="snp,scan"))
  
  unlink(file)
}

test_equalGenotypeDim <- function() {
  file <- tempfile()
  gds <- createfn.gds(file)
  snp <- 1:78
  chrom <- rep(1:26, each=3)
  pos <- rep(1001:1026, 3)
  a <- rep("A", 78)
  b <- rep("G", 78)
  alleles <- paste(a, b, sep="/")
  samp <- 1231:(1231+78-1)
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  add.gdsn(gds, "snp.chromosome", chrom)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  add.gdsn(gds, "genotype", t(geno), storage="bit2") # transpose of geno here for genotypeDim=scan,snp
  closefn.gds(gds)
  
  obj <- GdsGenotypeReader(file, genotypeDim="scan,snp")
  checkIdentical(obj@genotypeDim, "scan,snp")
  checkIdentical(snp, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(a, getAlleleA(obj))
  checkIdentical(b, getAlleleB(obj))
  checkIdentical(samp, getScanID(obj))
  geno[geno == 3] <- NA
  checkIdentical(geno, getGenotype(obj))
  checkIdentical(t(geno), getGenotype(obj, transpose=TRUE))
  
  sel <- samp %in% sample(samp, 3)
  checkIdentical(samp[sel], getScanID(obj, sel))
  close(obj)
  
  
  # this should raise an exception - snp and scan dimensions are equal
  checkException(GdsGenotypeReader(file))
  # this one will not raise an exception -- it would be the user's fault
  #checkException(GdsGenotypeReader(file, genotypeDim="snp,scan"))
  
  unlink(file)
  
}

test_indels <- function() {
  file <- tempfile()
  gds <- createfn.gds(file)
  snp <- 1:4
  chrom <- rep(1, 4)
  pos <- 101:104
  a <- c("A", "A", "AA", "AA")
  b <- c("G", "GG", "G", "GG")
  alleles <- paste(a, b, sep="/")
  samp <- 1231:1235
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  add.gdsn(gds, "snp.chromosome", chrom)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  add.gdsn(gds, "genotype", geno, storage="bit2")
  closefn.gds(gds)
  
  obj <- GdsGenotypeReader(file)
  checkIdentical(a, getAlleleA(obj))
  checkIdentical(b, getAlleleB(obj))
  close(obj)
  unlink(file)
}


test_logicalIndex <- function() {
  checkIdentical(c(rep(TRUE, 10), rep(FALSE, 10)),
                 GWASTools:::.logicalIndex(1:10, 20))
  checkException(GWASTools:::.logicalIndex(c(TRUE, FALSE), 3))
}

test_GdsGenotypeReader_index <- function() {  
  file <- tempfile()
  gds <- createfn.gds(file)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  a <- rep("A", 260)
  b <- rep("G", 260)
  alleles <- paste(a, b, sep="/")
  samp <- 1231:1235
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  add.gdsn(gds, "snp.chromosome", chrom)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  add.gdsn(gds, "genotype", geno, storage="bit2")
  closefn.gds(gds)
  
  obj <- GdsGenotypeReader(file)
  checkIdentical(snp[1:10], getSnpID(obj, index=1:10))
  checkIdentical(chrom[1:10], getChromosome(obj, index=1:10))
  checkIdentical(samp[1:2], getScanID(obj, index=c(rep(TRUE, 2), rep(FALSE, 3))))
  checkException(getScanID(obj, index=rep(TRUE, 6)))

  geno[geno == 3] <- NA
  checkIdentical(geno[1:10,1:2],
                 getGenotypeSelection(obj, snp=1:10, scan=1:2))
  checkIdentical(t(geno[1:10,1:2]),
                 getGenotypeSelection(obj, snp=1:10, scan=1:2, transpose=TRUE))
  checkIdentical(geno[1:10,1:2],
                 getGenotypeSelection(obj,
                                      snp=c(rep(TRUE, 10), rep(FALSE, 250)),
                                      scan=c(rep(TRUE, 2), rep(FALSE, 3))))
  
  close(obj)
  unlink(file)
}
