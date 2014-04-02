test_gdsSubsetCheck <- function(){
  
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
  # make sure it works  
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                 sample.include=sampID[sampsel],
                 snp.include=snpID[snpsel])

  # try a different sampsel
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                                sample.include=sampID[sampsel][-1],
                                snp.include=snpID[snpsel]))
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                                sample.include=sort(sample(sampID, length(sampsel))),
                                snp.include=snpID[snpsel]))
  
  # try a different snpsel
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                                sample.include=sampID[sampsel],
                                snp.include=snpID[snpsel][-1]))
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                                sample.include=sampID[sampsel],
                                snp.include=sort(sample(snpID, length(sampsel)))))
  
  unlink(gdsfile)
  unlink(subfile)
}


test_gdsSubsetCheck_geno <- function(){
  
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
  # make sure it works  
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                 sample.include=sampID[sampsel],
                 snp.include=snpID[snpsel])
  
  # change the sub file
  gds.sub <- openfn.gds(subfile, readonly=FALSE)
  geno.node <- index.gdsn(gds.sub, "genotype")
  i <- sample(10, 1)
  vals <- read.gdsn(geno.node, start=c(1, i), count=c(-1, 1))
  # change a missing value to a genotype
  j <- sample(which(vals == 3), 1)
  vals.new <- vals
  vals.new[j] <- 0
  write.gdsn(geno.node, vals.new, start=c(1, i), count=c(-1, 1))
  closefn.gds(gds.sub)

  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                 sample.include=sampID[sampsel],
                 snp.include=snpID[snpsel]))
  
  # put it back
  gds.sub <- openfn.gds(subfile, readonly=FALSE)
  geno.node <- index.gdsn(gds.sub, "genotype")
  write.gdsn(geno.node, vals, start=c(1, i), count=c(-1, 1))
  closefn.gds(gds.sub)
  
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                 sample.include=sampID[sampsel],
                 snp.include=snpID[snpsel])

  # change a genotype to a different genotype
  gds.sub <- openfn.gds(subfile, readonly=FALSE)
  geno.node <- index.gdsn(gds.sub, "genotype")
  i <- sample(10, 1)
  vals <- read.gdsn(geno.node, start=c(1, i), count=c(-1, 1))
  j <- sample(which(vals < 3), 1)
  vals.new <- vals
  vals.new[j] <- ifelse(vals[j] == 0, 2, vals[j]-1)
  write.gdsn(geno.node, vals.new, start=c(1, i), count=c(-1, 1))
  closefn.gds(gds.sub)
  
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                                sample.include=sampID[sampsel],
                                snp.include=snpID[snpsel]))
 
  unlink(gdsfile)
  unlink(subfile)
}


test_gdsSubsetCheck_attribute <- function(){
  
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
  # make sure it works  
  gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                 sample.include=sampID[sampsel],
                 snp.include=snpID[snpsel])
  

  # change an attribute in the sub file
  gds.sub <- openfn.gds(subfile, readonly=FALSE)
  geno.node <- index.gdsn(gds.sub, "genotype")
  put.attr.gdsn(geno.node, "test", 2)
  closefn.gds(gds.sub)
  
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                 sample.include=sampID[sampsel],
                 snp.include=snpID[snpsel]))
  
  # change the attribute in the parent file to a different value
  gds.parent <- openfn.gds(gdsfile, readonly=FALSE)
  geno.node <- index.gdsn(gds.parent, "genotype")
  put.attr.gdsn(geno.node, "test", 3)
  closefn.gds(gds.parent)
  
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile,
                                sample.include=sampID[sampsel],
                                snp.include=snpID[snpsel]))
  
  unlink(gdsfile)
  unlink(subfile)
}


test_gdsSubsetCheck_storage <- function(){
  
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
  geno <- matrix(runif(nsnp*nsamp, min=0, max=2),
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
  checkException(gdsSubsetCheck(parent.gds=gdsfile, sub.gds=subfile, snp.include=snpID[snpsel], sample.include=sampID[sampsel], sub.storage="bit2"))
  
  unlink(subfile)
  unlink(gdsfile)
  
}

