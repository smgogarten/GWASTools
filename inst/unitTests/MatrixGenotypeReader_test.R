test_MatrixGenotypeReader <- function() {
  snpID <- 1:100
  chrom <- rep(1:20, each=5)
  pos <- 1001:1100
  scanID <- 1:20
  geno <- matrix(sample(c(0,1,2,NA), 2000, replace=TRUE), nrow=100, ncol=20)

  mgr <- MatrixGenotypeReader(genotype=geno, snpID=snpID, chromosome=chrom, position=pos, scanID=scanID)

  checkEquals(snpID, getSnpID(mgr))
  checkEquals(scanID, getScanID(mgr))
  checkEquals(chrom, getChromosome(mgr))
  checkEquals(pos, getPosition(mgr))
  checkEquals(geno, getGenotype(mgr))
  checkEquals(geno[2:10, 5:6], getGenotype(mgr, snp=c(2,9), scan=c(5,2)))
  checkEquals(geno[2:10,], getGenotype(mgr, snp=c(2,9), scan=c(1,-1)))
  checkEquals(geno[2:10, 5:6], getGenotypeSelection(mgr, snp=2:10, scan=5:6, use.names=FALSE))
  checkEquals(snpID[1:5], getSnpID(mgr, index=1:5))

  checkEquals(geno[1,,drop=FALSE], getGenotype(mgr, snp=c(1,1), drop=FALSE))
  checkEquals(geno[1,,drop=FALSE], getGenotypeSelection(mgr, snp=1, drop=FALSE, use.names=FALSE))

  dimnames(geno) <- list(snpID, scanID)
  checkEquals(geno, getGenotype(mgr, use.names=TRUE))
  checkEquals(geno[2:10, 5:6], getGenotype(mgr, snp=c(2,9), scan=c(5,2), use.names=TRUE))
  checkEquals(geno[2:10, 5:6], getGenotypeSelection(mgr, snp=2:10, scan=5:6, use.names=TRUE))
  
  checkException(MatrixGenotypeReader(genotype=geno, snpID=snpID, chromosome=chrom, position=pos, scanID=scanID[1:10]))
}
