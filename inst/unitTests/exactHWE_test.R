
.hweGenoData <- function(nsnp=100, nsamp=50) {
    geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE), nrow=nsnp, ncol=nsamp)
    geno[1,] <- 0 ## make one monomorphic
    mgr <- MatrixGenotypeReader(geno, snpID=1:nsnp, scanID=1:nsamp,
                                chromosome=rep(c(1L,23L), each=nsnp/2),
                                position=1:nsnp)

    scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=1:nsamp,
      sex=sample(c("M","F"), nsamp, replace=TRUE), stringsAsFactors=FALSE))

    GenotypeData(mgr, scanAnnot=scanAnnot)
}


.checkHWE <- function(genoData, chromosome, scan.exclude=NULL) {
    chr <- range(which(getChromosome(genoData) == chromosome))
    hwe <- exactHWE(genoData, scan.exclude=scan.exclude, snpStart=chr[1], snpEnd=chr[2])
    mono <- hwe$MAF == 0
    checkTrue(all(is.na(hwe$pval[mono])))
    checkTrue(!all(is.na(hwe$pval)))

    cnt <- hwe[,c("nAA", "nAB", "nBB")]
    names(cnt) <- c("nAA", "nAa", "naa")
    tmp <- GWASExactHW::HWExact(cnt)
    checkEquals(hwe$pval[!mono], tmp[!mono])
}

test_HWE <- function() {
    genoData <- .hweGenoData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkHWE(genoData, chromosome=1, scan.exclude=scan.exclude)
    .checkHWE(genoData, chromosome=23, scan.exclude=scan.exclude)
    
    checkException(exactHWE(genoData)) # auto and X together
}


test_permuteGenotypes <- function() {
  
  genoData <- .hweGenoData()
  genotypes <- getGenotype(genoData)

  permuted <- GWASTools:::.permuteGenotypes(genotypes)

  # check same missing values
  checkEquals(is.na(permuted), is.na(genotypes))
  
  counts <- GWASTools:::.countGenotypes(genotypes)
  countsPermuted <- GWASTools:::.countGenotypes(permuted)
  
  checkEquals(2*counts[, "nAA"] + counts[, "nAa"],
              2*countsPermuted[, "nAA"] + countsPermuted[, "nAa"])  
  checkEquals(2*counts[, "naa"] + counts[, "nAa"],
              2*countsPermuted[, "naa"] + countsPermuted[, "nAa"])  
}

test_hwePermute <- function() {
  
  genoData <- .hweGenoData()
  
  snprange <- range(which(getChromosome(genoData) %in% 1:22))
  hwe <- exactHWE(genoData, snpStart=snprange[1], snpEnd=snprange[2])
  hwePermuted <- exactHWE(genoData, snpStart=snprange[1], snpEnd=snprange[2], permute=TRUE)

  checkEquals(hwe$MAF, hwePermuted$MAF)
  checkEquals(hwe$minor.allele, hwePermuted$minor.allele)
  checkEquals(is.na(hwe$f), is.na(hwePermuted$f))
  checkEquals(is.na(hwe$MAF), is.na(hwePermuted$MAF))
  
  snprange <- range(which(getChromosome(genoData) %in% 23))
  hwe <- exactHWE(genoData, snpStart=snprange[1], snpEnd=snprange[2])
  hwePermuted <- exactHWE(genoData, snpStart=snprange[1], snpEnd=snprange[2], permute=TRUE)
  
  checkEquals(hwe$MAF, hwePermuted$MAF)
  checkEquals(hwe$minor.allele, hwePermuted$minor.allele)
  checkEquals(is.na(hwe$f), is.na(hwePermuted$f))
  checkEquals(is.na(hwe$MAF), is.na(hwePermuted$MAF))
  
}
