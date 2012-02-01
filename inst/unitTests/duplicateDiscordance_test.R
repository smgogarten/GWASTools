test_duplicateDiscordance <- function() {
  # snp annotation
  snpID <- 1:10
  chrom <- rep(1L, 10)
  pos <- 101:110
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:6
  subjID <- c("a","b","c","b","b","a")
  scandf <- data.frame(scanID=scanID, subjID=subjID)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  
  # netCDF
  geno <- matrix(c(c(0,0,0,0,0,1,1,1,1,1),
                   c(1,1,1,1,1,2,2,2,2,2),
                   rep(0,10),
                   c(1,0,1,1,1,2,2,2,2,0),
                   c(1,1,0,1,0,2,2,2,2,0),
                   c(0,0,0,0,2,2,1,1,NA,1)), ncol=6)
  ncfile <- tempfile()
  ncdfCreate(snpdf, ncfile, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  # expected values
  a.exp <- c(0,0,0,0,1,1,0,0,0,0)
  b.exp <- c(0,2,2,0,2,0,0,0,0,2)
  subj.disc.exp <- c(0,1,1,0,2,1,0,0,0,1)
  tot.disc.exp <- c(0,2,2,0,3,1,0,0,0,2)
  npair.exp <- c(4,4,4,4,4,4,4,4,3,4)
  b.subj.exp <- matrix(c(0.0,0.2,0.3,0.2,0.0,0.3,0.3,0.3,0.0), ncol=3,
                       dimnames=list(c(2,4,5),c(2,4,5)))
  
  disc <- duplicateDiscordance(genoData, "subjID")
  checkIdentical(disc[[1]]$n.disc.subj, subj.disc.exp)
  checkIdentical(disc[[1]]$discordant, tot.disc.exp)
  checkIdentical(disc[[1]]$npair, npair.exp)
  checkEquals(disc[[2]]$b, b.subj.exp)

  # exclude some scans
  scan.exclude <- 5
  a.exp <- c(0,0,0,0,1,1,0,0,0,0)
  b.exp <- c(0,1,0,0,0,0,0,0,0,1)
  subj.disc.exp <- c(0,1,0,0,1,1,0,0,0,1)
  tot.disc.exp <- c(0,1,0,0,1,1,0,0,0,1)
  npair.exp <- c(2,2,2,2,2,2,2,2,1,2)
  b.subj.exp <- matrix(c(0.0,0.2,0.2,0.0), ncol=2,
                       dimnames=list(c(2,4),c(2,4)))
  
  disc <- duplicateDiscordance(genoData, "subjID", scan.exclude=scan.exclude)
  checkIdentical(disc[[1]]$n.disc.subj, subj.disc.exp)
  checkIdentical(disc[[1]]$discordant, tot.disc.exp)
  checkIdentical(disc[[1]]$npair, npair.exp)
  checkEquals(disc[[2]]$b, b.subj.exp)

  # exclude some snps
  snp.exclude <- c(2,10)
  a.exp <- c(0,0,0,1,1,0,0,0)
  b.exp <- c(0,2,0,2,0,0,0,0)
  subj.disc.exp <- c(0,1,0,2,1,0,0,0)
  tot.disc.exp <- c(0,2,0,3,1,0,0,0)
  npair.exp <- c(4,4,4,4,4,4,4,3)
  
  disc <- duplicateDiscordance(genoData, "subjID", snp.exclude=snp.exclude)
  checkIdentical(disc[[1]]$n.disc.subj, subj.disc.exp)
  checkIdentical(disc[[1]]$discordant, tot.disc.exp)
  checkIdentical(disc[[1]]$npair, npair.exp)

  # check that Y chrom SNPs for females are ignored
  snpAnnot$chromosome <- c(rep(1L, 2), rep(25L, 8))
  scanAnnot$sex <- "F"
  close(genoData)
  nc <- open.ncdf(ncfile, write=TRUE)
  put.var.ncdf(nc, "chromosome", snpAnnot$chromosome)
  close.ncdf(nc)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  snp.exclude <- 10
  a.exp <- c(0,0,0,0,0,0,0,0,0)
  b.exp <- c(0,2,0,0,0,0,0,0,0)
  subj.disc.exp <- c(0,1,0,0,0,0,0,0,0)
  tot.disc.exp <- c(0,2,0,0,0,0,0,0,0)
  npair.exp <- c(4,4,0,0,0,0,0,0,0)
  
  disc <- duplicateDiscordance(genoData, "subjID", snp.exclude=snp.exclude)
  checkIdentical(disc[[1]]$n.disc.subj, subj.disc.exp)
  checkIdentical(disc[[1]]$discordant, tot.disc.exp)
  checkIdentical(disc[[1]]$npair, npair.exp)
  
  close(genoData)
  file.remove(ncfile)
}
