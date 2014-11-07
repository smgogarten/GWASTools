test_gdsSubset <- function(){
  
  # simulated data
  gdsfile <- tempfile()
  subfile <- tempfile()
  
  gds <- createfn.gds(gdsfile)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  alleles <- c("A", "C", "T", "G")
  a <- sample(alleles, 260, replace=T)
  b <- sapply(a, function(x) sample(alleles[!(alleles %in% x)], 1), USE.NAMES=FALSE)
  alleles <- paste(a, b, sep="/")
  samp <- 1:20
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  node <- add.gdsn(gds, "snp.chromosome", chrom)
  put.attr.gdsn(node, "autosome.start", 1)
  put.attr.gdsn(node, "autosome.end", 22)
  put.attr.gdsn(node, "X", 23)
  put.attr.gdsn(node, "XY", 24)
  put.attr.gdsn(node, "Y", 25)
  put.attr.gdsn(node, "M", 26)
  put.attr.gdsn(node, "MT", 26)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  node <- add.gdsn(gds, "genotype", geno, storage="bit2")
  put.attr.gdsn(node, "snp.order")
  ## add in the same attributes later, when the function is working
  
  closefn.gds(gds)
  
  
  gds.orig <- GdsGenotypeReader(gdsfile)
  snpID <- getSnpID(gds.orig)
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- getScanID(gds.orig)
  sampsel <- sort(sample(1:length(sampID), 10))
  close(gds.orig)
  
  # subset
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel],
            snp.include=snpID[snpsel])
  
  
  # check attributes with gds
  gds1 <- openfn.gds(gdsfile)
  gds2 <- openfn.gds(subfile)
  checkTrue(setequal(ls.gdsn(gds1), ls.gdsn(gds2)))
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype")), get.attr.gdsn(index.gdsn(gds2, "genotype")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype"))$storage)
  closefn.gds(gds1)
  closefn.gds(gds2)
  # check the rest with GdsGentoypeReader
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID[snpsel])
  checkIdentical(getScanID(gds.sub), sampID[sampsel])
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, sampsel])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, sample.include=sampID[sampsel], snp.include=snpID[snpsel])
  unlink(subfile)
  
  
  # check with only sample include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel])
  
  # check
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID)
  checkIdentical(getScanID(gds.sub), sampID[sampsel])
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[, sampsel])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig))
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig))
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig))
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig))
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, sample.include=sampID[sampsel])
  unlink(subfile)
  
  # check with only snp include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            snp.include=snpID[snpsel])
  
  # check
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID[snpsel])
  checkIdentical(getScanID(gds.sub), sampID)
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, ])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel])
  unlink(subfile)  
  unlink(gdsfile)
}


test_gdsSubset_twoNodes <- function(){
  
  # simulated data
  gdsfile <- tempfile()
  subfile <- tempfile()
  
  gds <- createfn.gds(gdsfile)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  alleles <- c("A", "C", "T", "G")
  a <- sample(alleles, 260, replace=T)
  b <- sapply(a, function(x) sample(alleles[!(alleles %in% x)], 1), USE.NAMES=FALSE)
  alleles <- paste(a, b, sep="/")
  samp <- 1:20
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno1 <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                  nrow=nsnp, ncol=nsamp)
  geno2 <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                  nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  node <- add.gdsn(gds, "snp.chromosome", chrom)
  put.attr.gdsn(node, "autosome.start", 1)
  put.attr.gdsn(node, "autosome.end", 22)
  put.attr.gdsn(node, "X", 23)
  put.attr.gdsn(node, "XY", 24)
  put.attr.gdsn(node, "Y", 25)
  put.attr.gdsn(node, "M", 26)
  put.attr.gdsn(node, "MT", 26)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  node <- add.gdsn(gds, "genotype1", geno1, storage="bit2")
  put.attr.gdsn(node, "snp.order")
  node <- add.gdsn(gds, "genotype2", geno2, storage="bit2")
  put.attr.gdsn(node, "snp.order")
  
  
  closefn.gds(gds)
  
  #gds.orig <- GdsReader(gdsfile)
  #snpID <- getSnpID(gds.orig)
  #snpsel <- sort(sample(1:length(snpID), 50))
  #sampID <- getScanID(gds.orig)
  #sampsel <- sort(sample(1:length(sampID), 10))
  #close(gds.orig)
  
  snpID <- snp
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- samp
  sampsel <- sort(sample(1:length(sampID), 10))
  
  # subset
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel],
            snp.include=snpID[snpsel])
  
  
  # check node values with gdsfmt -- other attributes are checked are checked in other unit tests
  gds1 <- openfn.gds(gdsfile)
  gds2 <- openfn.gds(subfile)
  # genotype1  
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype1")), get.attr.gdsn(index.gdsn(gds2, "genotype1")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype1"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype1"))$storage)
  checkIdentical(read.gdsn(index.gdsn(gds1, "genotype1"))[snpsel, sampsel], read.gdsn(index.gdsn(gds2, "genotype1")))
  # genotype2
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype2")), get.attr.gdsn(index.gdsn(gds2, "genotype2")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype2"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype2"))$storage)
  checkIdentical(read.gdsn(index.gdsn(gds1, "genotype2"))[snpsel, sampsel], read.gdsn(index.gdsn(gds2, "genotype2")))
  # close files
  closefn.gds(gds1)
  closefn.gds(gds2)
  
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, sample.include=sampID[sampsel],
                 snp.include=snpID[snpsel])
  unlink(subfile)
  
  
  # check with only sample include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel])
  
  # check
  gds1 <- openfn.gds(gdsfile)
  gds2 <- openfn.gds(subfile)
  # genotype1  
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype1")), get.attr.gdsn(index.gdsn(gds2, "genotype1")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype1"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype1"))$storage)
  checkIdentical(read.gdsn(index.gdsn(gds1, "genotype1"))[, sampsel], read.gdsn(index.gdsn(gds2, "genotype1")))
  # genotype2
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype2")), get.attr.gdsn(index.gdsn(gds2, "genotype2")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype2"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype2"))$storage)
  checkIdentical(read.gdsn(index.gdsn(gds1, "genotype2"))[, sampsel], read.gdsn(index.gdsn(gds2, "genotype2")))
  # close files
  closefn.gds(gds1)
  closefn.gds(gds2)
  
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, sample.include=sampID[sampsel])
  unlink(subfile)
  
  
  # check with only snp include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            snp.include=snpID[snpsel])
  
  # check
  gds1 <- openfn.gds(gdsfile)
  gds2 <- openfn.gds(subfile)
  # genotype1  
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype1")), get.attr.gdsn(index.gdsn(gds2, "genotype1")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype1"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype1"))$storage)
  checkIdentical(read.gdsn(index.gdsn(gds1, "genotype1"))[snpsel, ], read.gdsn(index.gdsn(gds2, "genotype1")))
  # genotype2
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype2")), get.attr.gdsn(index.gdsn(gds2, "genotype2")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype2"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype2"))$storage)
  checkIdentical(read.gdsn(index.gdsn(gds1, "genotype2"))[snpsel, ], read.gdsn(index.gdsn(gds2, "genotype2")))
  # close files
  closefn.gds(gds1)
  closefn.gds(gds2)
  
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel])
  unlink(subfile)  
  unlink(gdsfile)
}

test_gdsSubset_scansnp <- function(){
  
  # simulated data
  gdsfile <- tempfile()
  subfile <- tempfile()
  
  gds <- createfn.gds(gdsfile)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  alleles <- c("A", "C", "T", "G")
  a <- sample(alleles, 260, replace=T)
  b <- sapply(a, function(x) sample(alleles[!(alleles %in% x)], 1), USE.NAMES=FALSE)
  alleles <- paste(a, b, sep="/")
  samp <- 21:40
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  node <- add.gdsn(gds, "snp.chromosome", chrom)
  put.attr.gdsn(node, "autosome.start", 1)
  put.attr.gdsn(node, "autosome.end", 22)
  put.attr.gdsn(node, "X", 23)
  put.attr.gdsn(node, "XY", 24)
  put.attr.gdsn(node, "Y", 25)
  put.attr.gdsn(node, "M", 26)
  put.attr.gdsn(node, "MT", 26)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  node <- add.gdsn(gds, "genotype", t(geno), storage="bit2")
  put.attr.gdsn(node, "sample.order")
  ## add in the same attributes later, when the function is working
  
  closefn.gds(gds)
  
  gds.orig <- GdsGenotypeReader(gdsfile)
  snpID <- getSnpID(gds.orig)
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- getScanID(gds.orig)
  sampsel <- sort(sample(1:length(sampID), 10))
  close(gds.orig)
  
  # subset
  for (block.size in c(5000, 20, 259, 1)) {
    subfile <- tempfile()
    gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
              sample.include=sampID[sampsel],
              snp.include=snpID[snpsel], block.size=block.size)
    
    # check attributes with gds
    gds1 <- openfn.gds(gdsfile)
    gds2 <- openfn.gds(subfile)
    checkTrue(setequal(ls.gdsn(gds1), ls.gdsn(gds2)))
    checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype")), get.attr.gdsn(index.gdsn(gds2, "genotype")))
    checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype"))$storage)
    closefn.gds(gds1)
    closefn.gds(gds2)
    # check the rest with GdsGentoypeReader
    gds.sub <- GdsGenotypeReader(subfile)
    gds.orig <- GdsGenotypeReader(gdsfile)
    checkIdentical(getSnpID(gds.sub), snpID[snpsel])
    checkIdentical(getScanID(gds.sub), sampID[sampsel])
    checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, sampsel])
    checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
    checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
    checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
    checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
    checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
    checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
    checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
    checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
    checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
    #gdsSubsetCheck()
    
    close(gds.orig)
    close(gds.sub)
    gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel], sample.include=sampID[sampsel])
    unlink(subfile)
  }
  
  # check with only sample include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel])
  
  # check
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID)
  checkIdentical(getScanID(gds.sub), sampID[sampsel])
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[, sampsel])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig))
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig))
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig))
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig))
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  #gdsSubsetCheck()
  #gdsSubsetCheck()
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, sample.include=sampID[sampsel])
  unlink(subfile)
  
  # check with only snp include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            snp.include=snpID[snpsel])
  
  # check
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID[snpsel])
  checkIdentical(getScanID(gds.sub), sampID)
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, ])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  #gdsSubsetCheck()
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel])
  unlink(subfile)  
  unlink(gdsfile)
}




test_gdsSubset_scansnp_snpid <- function(){
  
  # simulated data
  gdsfile <- tempfile()
  subfile <- tempfile()
  
  gds <- createfn.gds(gdsfile)
  snp <- (1:260) + 50
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  alleles <- c("A", "C", "T", "G")
  a <- sample(alleles, 260, replace=T)
  b <- sapply(a, function(x) sample(alleles[!(alleles %in% x)], 1), USE.NAMES=FALSE)
  alleles <- paste(a, b, sep="/")
  samp <- 21:40
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  node <- add.gdsn(gds, "snp.chromosome", chrom)
  put.attr.gdsn(node, "autosome.start", 1)
  put.attr.gdsn(node, "autosome.end", 22)
  put.attr.gdsn(node, "X", 23)
  put.attr.gdsn(node, "XY", 24)
  put.attr.gdsn(node, "Y", 25)
  put.attr.gdsn(node, "M", 26)
  put.attr.gdsn(node, "MT", 26)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  node <- add.gdsn(gds, "genotype", t(geno), storage="bit2")
  put.attr.gdsn(node, "sample.order")
  ## add in the same attributes later, when the function is working
  
  closefn.gds(gds)
  
  gds.orig <- GdsGenotypeReader(gdsfile)
  snpID <- getSnpID(gds.orig)
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- getScanID(gds.orig)
  sampsel <- sort(sample(1:length(sampID), 10))
  close(gds.orig)
  
  # subset
  
    subfile <- tempfile()
    gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
              sample.include=sampID[sampsel],
              snp.include=snpID[snpsel], block.size=25)
    
    # check attributes with gds
    gds1 <- openfn.gds(gdsfile)
    gds2 <- openfn.gds(subfile)
    checkTrue(setequal(ls.gdsn(gds1), ls.gdsn(gds2)))
    checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype")), get.attr.gdsn(index.gdsn(gds2, "genotype")))
    checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype"))$storage)
    closefn.gds(gds1)
    closefn.gds(gds2)
    # check the rest with GdsGentoypeReader
    gds.sub <- GdsGenotypeReader(subfile)
    gds.orig <- GdsGenotypeReader(gdsfile)
    checkIdentical(getSnpID(gds.sub), snpID[snpsel])
    checkIdentical(getScanID(gds.sub), sampID[sampsel])
    checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, sampsel])
    checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
    checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
    checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
    checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
    checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
    checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
    checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
    checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
    checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
    #gdsSubsetCheck()
    
    close(gds.orig)
    close(gds.sub)
    gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel], sample.include=sampID[sampsel])
    unlink(subfile)
  
  # check with only sample include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel])
  
  # check
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID)
  checkIdentical(getScanID(gds.sub), sampID[sampsel])
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[, sampsel])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig))
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig))
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig))
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig))
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  #gdsSubsetCheck()
  #gdsSubsetCheck()
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, sample.include=sampID[sampsel])
  unlink(subfile)
  
  # check with only snp include
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            snp.include=snpID[snpsel])
  
  # check
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID[snpsel])
  checkIdentical(getScanID(gds.sub), sampID)
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, ])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  #gdsSubsetCheck()
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel])
  unlink(subfile)  
  unlink(gdsfile)
}


test_gdsSubset_storage <- function(){
  
  # simulated data
  gdsfile <- tempfile()
  subfile <- tempfile()
  
  gds <- createfn.gds(gdsfile)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  alleles <- c("A", "C", "T", "G")
  a <- sample(alleles, 260, replace=T)
  b <- sapply(a, function(x) sample(alleles[!(alleles %in% x)], 1), USE.NAMES=FALSE)
  alleles <- paste(a, b, sep="/")
  samp <- 21:40
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(-1:2, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  node <- add.gdsn(gds, "snp.chromosome", chrom)
  put.attr.gdsn(node, "autosome.start", 1)
  put.attr.gdsn(node, "autosome.end", 22)
  put.attr.gdsn(node, "X", 23)
  put.attr.gdsn(node, "XY", 24)
  put.attr.gdsn(node, "Y", 25)
  put.attr.gdsn(node, "M", 26)
  put.attr.gdsn(node, "MT", 26)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  node <- add.gdsn(gds, "genotype", geno, storage="float")
  put.attr.gdsn(node, "snp.order")
  put.attr.gdsn(node, "missing.value", -1)
  ## add in the same attributes later, when the function is working
  
  closefn.gds(gds)
  
  gds.orig <- GdsGenotypeReader(gdsfile)
  snpID <- getSnpID(gds.orig)
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- getScanID(gds.orig)
  sampsel <- sort(sample(1:length(sampID), 10))
  close(gds.orig)
  
  # subset
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel],
            snp.include=snpID[snpsel], sub.storage="bit2")
  
  # check attributes with gds
  gds1 <- openfn.gds(gdsfile)
  gds2 <- openfn.gds(subfile)
  checkTrue(setequal(ls.gdsn(gds1), ls.gdsn(gds2)))
  #checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype")), get.attr.gdsn(index.gdsn(gds2, "genotype")))
  checkEquals(get.attr.gdsn(index.gdsn(gds2, "genotype"))[["missing.value"]], 3)
  checkIdentical(objdesp.gdsn(index.gdsn(gds2, "genotype"))$storage, "Bit2")
  closefn.gds(gds1)
  closefn.gds(gds2)
  # check the rest with GdsGentoypeReader
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID[snpsel])
  checkIdentical(getScanID(gds.sub), sampID[sampsel])
  checkEquals(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, sampsel]) # not identical, different storage types
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  #gdsSubsetCheck()
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel], sample.include=sampID[sampsel], sub.storage="bit2")
  unlink(subfile)
  unlink(gdsfile)
  
}

test_gdsSubset_attribute <- function(){
  
  # simulated data
  gdsfile <- tempfile()
  subfile <- tempfile()
  
  gds <- createfn.gds(gdsfile)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  alleles <- c("A", "C", "T", "G")
  a <- sample(alleles, 260, replace=T)
  b <- sapply(a, function(x) sample(alleles[!(alleles %in% x)], 1), USE.NAMES=FALSE)
  alleles <- paste(a, b, sep="/")
  samp <- 21:40
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  node <- add.gdsn(gds, "snp.chromosome", chrom)
  put.attr.gdsn(node, "autosome.start", 1)
  put.attr.gdsn(node, "autosome.end", 22)
  put.attr.gdsn(node, "X", 23)
  put.attr.gdsn(node, "XY", 24)
  put.attr.gdsn(node, "Y", 25)
  put.attr.gdsn(node, "M", 26)
  put.attr.gdsn(node, "MT", 26)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  node <- add.gdsn(gds, "genotype", geno, storage="bit2")
  put.attr.gdsn(node, "snp.order")
  add.gdsn(gds, "test", 5)
  
  closefn.gds(gds)
  
  gds.orig <- GdsGenotypeReader(gdsfile)
  snpID <- getSnpID(gds.orig)
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- getScanID(gds.orig)
  sampsel <- sort(sample(1:length(sampID), 10))
  close(gds.orig)
  
  # subset
  subfile <- tempfile()
  gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel],
            snp.include=snpID[snpsel])
  
  
  # check attributes with gds
  gds1 <- openfn.gds(gdsfile)
  gds2 <- openfn.gds(subfile)
  checkTrue(setequal(ls.gdsn(gds1), c(ls.gdsn(gds2), "test")))
  checkIdentical(get.attr.gdsn(index.gdsn(gds1, "genotype")), get.attr.gdsn(index.gdsn(gds2, "genotype")))
  checkIdentical(objdesp.gdsn(index.gdsn(gds1, "genotype"))$storage, objdesp.gdsn(index.gdsn(gds2, "genotype"))$storage)
  closefn.gds(gds1)
  closefn.gds(gds2)
  # check the rest with GdsGentoypeReader
  gds.sub <- GdsGenotypeReader(subfile)
  gds.orig <- GdsGenotypeReader(gdsfile)
  checkIdentical(getSnpID(gds.sub), snpID[snpsel])
  checkIdentical(getScanID(gds.sub), sampID[sampsel])
  checkIdentical(getGenotype(gds.sub), getGenotype(gds.orig)[snpsel, sampsel])
  checkIdentical(getChromosome(gds.sub), getChromosome(gds.orig)[snpsel])
  checkIdentical(getPosition(gds.sub), getPosition(gds.orig)[snpsel])
  checkIdentical(getAlleleA(gds.sub), getAlleleA(gds.orig)[snpsel])
  checkIdentical(getAlleleB(gds.sub), getAlleleB(gds.orig)[snpsel])
  checkIdentical(autosomeCode(gds.sub), autosomeCode(gds.orig))
  checkIdentical(XchromCode(gds.sub), XchromCode(gds.orig))
  checkIdentical(YchromCode(gds.sub), YchromCode(gds.orig))
  checkIdentical(XYchromCode(gds.sub), XYchromCode(gds.orig))
  checkIdentical(MchromCode(gds.sub), MchromCode(gds.orig))
  
  close(gds.orig)
  close(gds.sub)
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, sample.include=sampID[sampsel], snp.include=snpID[snpsel])
  unlink(subfile)
  
  unlink(gdsfile)
}


test_gdsSubset_noOrder <- function(){
  
  # simulated data
  gdsfile <- tempfile()
  subfile <- tempfile()
  
  gds <- createfn.gds(gdsfile)
  snp <- 1:260
  chrom <- rep(1:26, each=10)
  pos <- rep(1001:1026, 10)
  alleles <- c("A", "C", "T", "G")
  a <- sample(alleles, 260, replace=T)
  b <- sapply(a, function(x) sample(alleles[!(alleles %in% x)], 1), USE.NAMES=FALSE)
  alleles <- paste(a, b, sep="/")
  samp <- 1:20
  nsnp <- length(snp)
  nsamp <- length(samp)
  geno <- matrix(sample(0:3, nsnp*nsamp, replace=TRUE),
                 nrow=nsnp, ncol=nsamp)
  add.gdsn(gds, "snp.id", snp)
  node <- add.gdsn(gds, "snp.chromosome", chrom)
  put.attr.gdsn(node, "autosome.start", 1)
  put.attr.gdsn(node, "autosome.end", 22)
  put.attr.gdsn(node, "X", 23)
  put.attr.gdsn(node, "XY", 24)
  put.attr.gdsn(node, "Y", 25)
  put.attr.gdsn(node, "M", 26)
  put.attr.gdsn(node, "MT", 26)
  add.gdsn(gds, "snp.position", pos)
  add.gdsn(gds, "snp.allele", alleles)
  add.gdsn(gds, "sample.id", samp)
  node <- add.gdsn(gds, "genotype", geno, storage="bit2")
  ## add in the same attributes later, when the function is working
  
  closefn.gds(gds)
  
  
  gds.orig <- GdsGenotypeReader(gdsfile)
  snpID <- getSnpID(gds.orig)
  snpsel <- sort(sample(1:length(snpID), 50))
  sampID <- getScanID(gds.orig)
  sampsel <- sort(sample(1:length(sampID), 10))
  close(gds.orig)
  
  # subset
  subfile <- tempfile()
  checkException(gdsSubset(parent.gds=gdsfile, sub.gds=subfile,
            sample.include=sampID[sampsel],
            snp.include=snpID[snpsel]))
}
