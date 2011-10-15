test_discordantPair <- function() {
  # first set
  # snp annotation
  snpID <- 1:10
  chrom <- rep(1L, 10)
  pos <- 101:110
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:5
  subjID <- c("a","b","c","d","e")
  scandf <- data.frame(scanID=scanID, subjID=subjID)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  
  # netCDF
  geno <- matrix(c(rep(0,5), rep(1,5),
                   rep(1,5), rep(2,5),
                   rep(2,5), rep(0,5),
                   rep(0,10),
                   rep(0,10)), ncol=5)
  ncfile1 <- tempfile()
  ncdfCreate(snpdf, ncfile1, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile1, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc1 <- NcdfGenotypeReader(ncfile1)
  genoData1 <- GenotypeData(nc1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  # second set
  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- as.integer(c(101,103,105,107,109))
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:4
  subjID <- c("c","f","b","a")
  scandf <- data.frame(scanID=scanID, subjID=subjID)
  scanAnnot <- ScanAnnotationDataFrame(scandf)

  # netCDF
  geno <- matrix(c(2,1,NA,1,0,
                   rep(0,5),
                   1,1,0,2,2,
                   2,0,0,1,1), ncol=4)
  ncfile2 <- tempfile()
  ncdfCreate(snpdf, ncfile2, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile2, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc2 <- NcdfGenotypeReader(ncfile2)
  genoData2 <- GenotypeData(nc2, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  # expected output (TRUE for discordance)
  # sample "a"
  a.exp <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  a.nm <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  # sample "b"
  b.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  b.nm <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  # sample "c"
  c.exp <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
  c.nm <- c(TRUE, TRUE, FALSE, TRUE, TRUE)

  # id matching
  a.scanID1 <- 1
  a.scanID2 <- 4
  b.scanID1 <- 2
  b.scanID2 <- 3
  c.scanID1 <- 3
  c.scanID2 <- 1
  snpID1 <- c(1,3,5,7,9)
  snpID2 <- 1:5

  # test
  a <- GWASTools:::discordantPair(genoData1, a.scanID1, snpID1,
                       genoData2, a.scanID2, snpID2)
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- GWASTools:::discordantPair(genoData1, b.scanID1, snpID1,
                       genoData2, b.scanID2, snpID2)
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
  c <- GWASTools:::discordantPair(genoData1, c.scanID1, snpID1,
                       genoData2, c.scanID2, snpID2)
  checkIdentical(c$discordant, c.exp)
  checkIdentical(c$nonmissing, c.nm)

  # what if we rearrange the order of SNPs?
  snpID1 <- c(3,5,1,9,7)
  snpID2 <- c(2,3,1,5,4)
  a.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  a.nm <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  b.exp <- c(FALSE, TRUE, FALSE, FALSE, FALSE)
  b.nm <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  c.exp <- c(TRUE, FALSE, FALSE, FALSE, TRUE)
  c.nm <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
  a <- GWASTools:::discordantPair(genoData1, a.scanID1, snpID1,
                       genoData2, a.scanID2, snpID2)
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- GWASTools:::discordantPair(genoData1, b.scanID1, snpID1,
                       genoData2, b.scanID2, snpID2)
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
  c <- GWASTools:::discordantPair(genoData1, c.scanID1, snpID1,
                       genoData2, c.scanID2, snpID2)
  checkIdentical(c$discordant, c.exp)
  checkIdentical(c$nonmissing, c.nm)

  # check error conditions
  # snpID1 is too long
  checkException(GWASTools:::discordantPair(genoData1, a.scanID1, 1:10,
                                 genoData2, a.scanID2, snpID2))
  # snpID1 has wrong values
  checkException(GWASTools:::discordantPair(genoData1, a.scanID1, 16:20,
                                 genoData2, a.scanID2, snpID2))
  # scanID1 has wrong value
  checkException(GWASTools:::discordantPair(genoData1, 10, snpID1,
                                 genoData2, a.scanID2, snpID2))
  
  file.remove(ncfile1, ncfile2)
}



test_duplicateDiscordanceAcrossDatasets <- function() {
  # test 1: only one scan per subject in both datasets
  # first set
  # snp annotation
  snpID <- 1:10
  chrom <- rep(1L, 10)
  pos <- 101:110
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:5
  subjID <- c("a","b","c","d","e")
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  
  # netCDF
  geno <- matrix(c(rep(0,5), rep(1,5),
                   rep(1,5), rep(2,5),
                   rep(2,5), rep(0,5),
                   rep(0,10),
                   rep(0,10)), ncol=5)
  ncfile1 <- tempfile()
  ncdfCreate(snpdf, ncfile1, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile1, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc1 <- NcdfGenotypeReader(ncfile1)
  genoData1 <- GenotypeData(nc1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  # second set
  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- as.integer(c(101,103,105,107,109))
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:5
  subjID <- c("c","f","b","a",NA)
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot <- ScanAnnotationDataFrame(scandf)

  # netCDF
  geno <- matrix(c(2,NA,1,1,0,
                   rep(0,5),
                   1,1,0,2,NA,
                   2,0,0,1,NA,
                   rep(0,5)), ncol=5)
  ncfile2 <- tempfile()
  ncdfCreate(snpdf, ncfile2, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile2, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc2 <- NcdfGenotypeReader(ncfile2)
  genoData2 <- GenotypeData(nc2, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  # expected output (TRUE for discordance)
  a.exp <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  b.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  c.exp <- c(FALSE, FALSE, TRUE, TRUE, FALSE)
  snp.exp <- data.frame(discordant=c(1,0,2,1,0), npair=c(3,2,3,3,1),
    n.disc.subj=c(1,0,2,1,0), discord.rate=c(1,0,2,1,0)/c(3,2,3,3,1))
  row.names(snp.exp) <- c("rs101","rs103","rs105","rs107","rs109")
  subj.exp <- list(a=matrix(1/4, 1, 1, dimnames=list(1,4)),
                   b=matrix(1/4, 1, 1, dimnames=list(2,3)),
                   c=matrix(2/4, 1, 1, dimnames=list(3,1)))
  
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2))
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)

  
  # test 2: add a second scan for subject "c" to dataset 2
  # scan annotation
  scanID <- 1:5
  subjID <- c("c","f","b","a","c")
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot <- ScanAnnotationDataFrame(scandf)

  # netCDF
  geno <- matrix(c(2,1,1,1,0,
                   rep(0,5),
                   1,1,0,2,2,
                   2,0,0,1,1,
                   2,2,2,1,0), ncol=5)
  ncdfCreate(snpdf, ncfile2, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile2, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc2 <- NcdfGenotypeReader(ncfile2)
  genoData2 <- GenotypeData(nc2, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  # expected output (TRUE for discordance)
  a.exp <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  b.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  c1.exp <- c(FALSE, TRUE, TRUE, TRUE, FALSE)
  c2.exp <- c(FALSE, FALSE, FALSE, TRUE, FALSE)
  snp.exp <- data.frame(discordant=c(1,1,2,2,0), npair=c(4,4,4,4,4),
    n.disc.subj=c(1,1,2,1,0), discord.rate=c(1,1,2,2,0)/c(4,4,4,4,4))
  row.names(snp.exp) <- c("rs101","rs103","rs105","rs107","rs109")
  subj.exp <- list(a=matrix(1/5, 1, 1, dimnames=list(1,4)),
                   b=matrix(1/5, 1, 1, dimnames=list(2,3)),
                   c=matrix(c(3/5, 1/5), 1, 2, dimnames=list(3,c(1,5))))

  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2))
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)

  
  # select particular snps to include
  snp.include <- c("rs101","rs103","rs105")
  snp.exp <- data.frame(discordant=c(1,1,2), npair=c(4,4,4),
    n.disc.subj=c(1,1,2), discord.rate=c(1,1,2)/c(4,4,4))
  row.names(snp.exp) <- c("rs101","rs103","rs105")
  subj.exp <- list(a=matrix(1/3, 1, 1, dimnames=list(1,4)),
                   b=matrix(1/3, 1, 1, dimnames=list(2,3)),
                   c=matrix(c(2/3, 0), 1, 2, dimnames=list(3,c(1,5))))
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2), snp.include=snp.include)
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)

  
  # change the order of the snps
  snp.include <- c("rs103","rs105","rs101")
  snp.exp <- data.frame(discordant=c(1,2,1), npair=c(4,4,4),
    n.disc.subj=c(1,2,1), discord.rate=c(1,2,1)/c(4,4,4))
  row.names(snp.exp) <- c("rs103","rs105","rs101")
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2), snp.include=snp.include)
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)

  
  # exclude scans
  excl1 <- 2
  excl2 <- 5
  snp.exp <- data.frame(discordant=c(1,1,1,1,0), npair=c(2,2,2,2,2),
    n.disc.subj=c(1,1,1,1,0), discord.rate=c(1,1,1,1,0)/c(2,2,2,2,2))
  row.names(snp.exp) <- c("rs101","rs103","rs105","rs107","rs109")
  subj.exp <- list(a=matrix(1/5, 1, 1, dimnames=list(1,4)),
                   c=matrix(3/5, 1, 1, dimnames=list(3,1)))

  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2), scan.exclude1=excl1, scan.exclude2=excl2)
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)
  
  
  # error check - snps not in dataset 2
  snp.include <- c("rs102", "rs104")
  checkException({
    discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
      rep("subjID", 2), rep("rsID", 2), snp.include=snp.include)
  })

  
  # error check - datasets with no common snps (expect NULL)
  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- 111:115
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:4
  subjID <- c("c","f","b","a")
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot <- ScanAnnotationDataFrame(scandf)

  # netCDF
  geno <- matrix(c(rep(0,5), rep(0,5), rep(0,5), rep(0,5)), ncol=4)
  ncdfCreate(snpdf, ncfile2, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile2, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc2 <- NcdfGenotypeReader(ncfile2)
  genoData2 <- GenotypeData(nc2, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
      rep("subjID", 2), rep("rsID", 2))
  checkIdentical(discord, NULL)

  
  # error check - datasets with no common subjects (expect NULL)
  # snp annotation
  snpID <- 1:10
  chrom <- rep(1L, 10)
  pos <- 101:110
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:5
  subjID <- c("a","b","c","d",NA)
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  
  # netCDF
  geno <- matrix(c(rep(0,5), rep(1,5),
                   rep(1,5), rep(2,5),
                   rep(2,5), rep(0,5),
                   rep(0,10),
                   rep(0,10)), ncol=5)
  ncdfCreate(snpdf, ncfile1, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile1, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc1 <- NcdfGenotypeReader(ncfile1)
  genoData1 <- GenotypeData(nc1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- as.integer(c(101,103,105,107,109))
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:4
  subjID <- c("g","f",NA,NA)
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot <- ScanAnnotationDataFrame(scandf)

  # netCDF
  geno <- matrix(c(rep(0,5), rep(0,5), rep(0,5), rep(0,5)), ncol=4)
  ncdfCreate(snpdf, ncfile2, n.samples=nrow(scandf))
  nc <- open.ncdf(ncfile2, write=TRUE)
  put.var.ncdf(nc, "sampleID", scanID)
  put.var.ncdf(nc, "genotype", geno)
  close.ncdf(nc)
  nc2 <- NcdfGenotypeReader(ncfile2)
  genoData2 <- GenotypeData(nc2, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
  
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
      rep("subjID", 2), rep("rsID", 2))
  checkIdentical(discord, NULL)

  close(genoData1)
  close(genoData2)
  file.remove(ncfile1, ncfile2)
}
