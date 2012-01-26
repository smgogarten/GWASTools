test_ncdfSubset <- function() {
  # simulated data 
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)
  snpID <- getSnpID(nc)
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- getScanID(nc)
  sampsel <- sort(sample(1:length(sampID), 10))
  geno <- getGenotype(nc)
  chrom <- getChromosome(nc)
  pos <- getPosition(nc)
  close(nc)

  # check subset
  subfile <- tempfile()
  ncdfSubset(parent.ncdf=ncfile, sub.ncdf=subfile,
             sample.include=sampID[sampsel], snp.include=snpID[snpsel])
  subnc <- NcdfGenotypeReader(subfile)
  subsnp <- getSnpID(subnc)
  checkIdentical(subsnp, snpID[snpsel])
  subsamp <- getScanID(subnc)
  checkIdentical(subsamp, sampID[sampsel])
  subgeno <- getGenotype(subnc)
  checkIdentical(subgeno, geno[snpsel, sampsel])
  subchrom <- getChromosome(subnc)
  checkIdentical(subchrom, chrom[snpsel])
  subpos <- getPosition(subnc)
  checkIdentical(subpos, pos[snpsel])
  ncdfSubsetCheck(parent.ncdf=ncfile, sub.ncdf=subfile,
                  sample.include=sampID[sampsel], snp.include=snpID[snpsel])
  close(subnc)
  file.remove(subfile)

  # check with only sample.include
  subfile <- tempfile()
  ncdfSubset(parent.ncdf=ncfile, sub.ncdf=subfile,
             sample.include=sampID[sampsel])
  subnc <- NcdfGenotypeReader(subfile)
  subsnp <- getSnpID(subnc)
  checkIdentical(subsnp, snpID)
  subsamp <- getScanID(subnc)
  checkIdentical(subsamp, sampID[sampsel])
  subgeno <- getGenotype(subnc)
  checkIdentical(subgeno, geno[, sampsel])
  subchrom <- getChromosome(subnc)
  checkIdentical(subchrom, chrom)
  subpos <- getPosition(subnc)
  checkIdentical(subpos, pos)
  ncdfSubsetCheck(parent.ncdf=ncfile, sub.ncdf=subfile,
             sample.include=sampID[sampsel])
  close(subnc)
  file.remove(subfile)
  
  # check with only snp.include
  subfile <- tempfile()
  ncdfSubset(parent.ncdf=ncfile, sub.ncdf=subfile,
             snp.include=snpID[snpsel])
  subnc <- NcdfGenotypeReader(subfile)
  subsnp <- getSnpID(subnc)
  checkIdentical(subsnp, snpID[snpsel])
  subsamp <- getScanID(subnc)
  checkIdentical(subsamp, sampID)
  subgeno <- getGenotype(subnc)
  checkIdentical(subgeno, geno[snpsel,])
  subchrom <- getChromosome(subnc)
  checkIdentical(subchrom, chrom[snpsel])
  subpos <- getPosition(subnc)
  checkIdentical(subpos, pos[snpsel])
  ncdfSubsetCheck(parent.ncdf=ncfile, sub.ncdf=subfile,
             snp.include=snpID[snpsel])
  close(subnc)
  file.remove(subfile)
  
  # check with both snp and sample include=NULL
  subfile <- tempfile()
  ncdfSubset(parent.ncdf=ncfile, sub.ncdf=subfile)
  subnc <- NcdfGenotypeReader(subfile)
  subsnp <- getSnpID(subnc)
  checkIdentical(subsnp, snpID)
  subsamp <- getScanID(subnc)
  checkIdentical(subsamp, sampID)
  subgeno <- getGenotype(subnc)
  checkIdentical(subgeno, geno)
  subchrom <- getChromosome(subnc)
  checkIdentical(subchrom, chrom)
  subpos <- getPosition(subnc)
  checkIdentical(subpos, pos)
  ncdfSubsetCheck(parent.ncdf=ncfile, sub.ncdf=subfile)
  close(subnc)
  file.remove(subfile)

  file.remove(ncfile)
}
