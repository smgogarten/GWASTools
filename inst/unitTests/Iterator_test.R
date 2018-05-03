
.testfile <- function() {
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
  file
}

.testannot <- function(x) {
    SnpAnnotationDataFrame(data.frame(
        snpID=getSnpID(x),
        chromosome=getChromosome(x),
        position=getPosition(x)
    ))
}

test_iterator <- function() {
  file <- .testfile()

  gdsobj <- GdsGenotypeReader(file)
  obj <- GenotypeData(gdsobj)
  it <- GenotypeBlockIterator(obj, snpBlock=100)
  checkEquals(1:100, currentFilter(it))
  checkEquals(1:100, getSnpID(it))
  checkEquals(getChromosome(obj)[1:100], getChromosome(it))
  checkEquals(getPosition(obj)[1:100], getPosition(it))
  checkEquals(getAlleleA(obj)[1:100], getAlleleA(it))
  checkEquals(getAlleleB(obj)[1:100], getAlleleB(it))
  checkEquals(getGenotypeSelection(obj)[1:100,], getGenotypeSelection(it))
  checkEquals(getGenotypeSelection(obj)[1:100,1:2], getGenotypeSelection(it, scan=1:2))
  checkTrue(iterateFilter(it))
  checkEquals(101:200, currentFilter(it))
  checkEquals(101:200, getSnpID(it))
  checkEquals(getGenotypeSelection(obj)[101:200,], getGenotypeSelection(it))
  checkTrue(iterateFilter(it))
  checkEquals(201:260, currentFilter(it))
  checkEquals(201:260, getSnpID(it))
  checkEquals(getGenotypeSelection(obj)[201:260,], getGenotypeSelection(it))
  checkTrue(!iterateFilter(it))

  close(obj)
  unlink(file)
}

test_iterator_snpAnnot <- function() {
  file <- .testfile()

  gdsobj <- GdsGenotypeReader(file)
  obj <- GenotypeData(gdsobj, snpAnnot=.testannot(gdsobj))
  it <- GenotypeBlockIterator(obj, snpBlock=100)
  checkEquals(1:100, getSnpVariable(it, "snpID"))
  checkTrue(iterateFilter(it))
  checkEquals(101:200, getSnpVariable(it, "snpID"))
  
  close(obj)
  unlink(file)
}

test_iterator_largeblock <- function() {
  file <- .testfile()

  gdsobj <- GdsGenotypeReader(file)
  obj <- GenotypeData(gdsobj)
  it <- GenotypeBlockIterator(obj, snpBlock=1000)
  checkEquals(1:260, currentFilter(it))
  checkEquals(1:260, getSnpID(it))
  checkTrue(!iterateFilter(it))
  
  close(obj)
  unlink(file)
}
