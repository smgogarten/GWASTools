
dp.test <- function(genoData1, scanID1, snpID1,
                    genoData2, scanID2, snpID2, ...) {
  geno1 <- GWASTools:::.selectGenotype(genoData1, scanID1, snpID1)
  geno2 <- GWASTools:::.selectGenotype(genoData2, scanID2, snpID2)
  GWASTools:::.discordantPair(geno1, geno2, ...)
}
  

test_discordantPair <- function() {
  # first set
  # snp annotation
  snpID <- 1:10
  chrom <- rep(1L, 10)
  pos <- 101:110
  alleleA <- rep("A",10)
  alleleB <- rep("G",10)
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      alleleA=alleleA, alleleB=alleleB,
                      stringsAsFactors=FALSE)
  snpAnnot1 <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:5
  subjID <- c("a","b","c","d","e")
  scandf <- data.frame(scanID=scanID, subjID=subjID)
  scanAnnot1 <- ScanAnnotationDataFrame(scandf)
  
  geno1 <- matrix(c(rep(0,5), rep(1,5),
                   rep(1,5), rep(2,5),
                   rep(2,5), rep(0,5),
                   rep(0,10),
                   rep(0,10)), ncol=5)
  mgr <- MatrixGenotypeReader(genotype=geno1, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData1 <- GenotypeData(mgr, snpAnnot=snpAnnot1, scanAnnot=scanAnnot1)

  # second set
  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- as.integer(c(101,103,105,107,109))
  alleleA <- rep("A",5)
  alleleB <- rep("G",5)
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      alleleA=alleleA, alleleB=alleleB,
                      stringsAsFactors=FALSE)
  snpAnnot2 <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:4
  subjID <- c("c","f","b","a")
  scandf <- data.frame(scanID=scanID, subjID=subjID)
  scanAnnot2 <- ScanAnnotationDataFrame(scandf)

  geno2 <- matrix(c(2,1,NA,1,0,
                   rep(0,5),
                   1,1,0,2,2,
                   2,0,0,1,1), ncol=4)
  mgr <- MatrixGenotypeReader(genotype=geno2, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData2 <- GenotypeData(mgr, snpAnnot=snpAnnot2, scanAnnot=scanAnnot2)

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
  a <- dp.test(genoData1, a.scanID1, snpID1,
               genoData2, a.scanID2, snpID2)
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- dp.test(genoData1, b.scanID1, snpID1,
               genoData2, b.scanID2, snpID2)
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
  c <- dp.test(genoData1, c.scanID1, snpID1,
               genoData2, c.scanID2, snpID2)
  checkIdentical(c$discordant, c.exp)
  checkIdentical(c$nonmissing, c.nm)

  # expected output for minor.allele.only
  major.genotype <- c("A/A","G/G","G/G","A/A","A/A")
  a.exp <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  a.nm <- c(TRUE, FALSE, FALSE, TRUE, TRUE)
  b.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  b.nm <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
  c.exp <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
  c.nm <- c(FALSE, TRUE, FALSE, TRUE, TRUE)
  a <- dp.test(genoData1, a.scanID1, snpID1,
               genoData2, a.scanID2, snpID2, major.genotype=major.genotype)
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- dp.test(genoData1, b.scanID1, snpID1,
               genoData2, b.scanID2, snpID2, major.genotype=major.genotype)
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
  c <- dp.test(genoData1, c.scanID1, snpID1,
               genoData2, c.scanID2, snpID2, major.genotype=major.genotype)
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
  a <- dp.test(genoData1, a.scanID1, snpID1,
               genoData2, a.scanID2, snpID2)
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- dp.test(genoData1, b.scanID1, snpID1,
               genoData2, b.scanID2, snpID2)
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
  c <- dp.test(genoData1, c.scanID1, snpID1,
               genoData2, c.scanID2, snpID2)
  checkIdentical(c$discordant, c.exp)
  checkIdentical(c$nonmissing, c.nm)

  # check error conditions
  # snpID1 is too long
  checkException(dp.test(genoData1, a.scanID1, 1:10,
                         genoData2, a.scanID2, snpID2))
  # snpID1 has wrong values
  checkException(dp.test(genoData1, a.scanID1, 16:20,
                         genoData2, a.scanID2, snpID2))
  # scanID1 has wrong value
  checkException(dp.test(genoData1, 10, snpID1,
                         genoData2, a.scanID2, snpID2))
  
  # check that Y chrom SNPs for females are ignored
  chrom <- c(rep(1L, 2), rep(25L, 8))
  snpAnnot1$chromosome <- chrom
  scanAnnot1$sex <- "F"
  mgr <- MatrixGenotypeReader(genotype=geno1, snpID=snpAnnot1$snpID,
    chromosome=chrom, position=snpAnnot1$pos, scanID=scanAnnot1$scanID)
  genoData1 <- GenotypeData(mgr, snpAnnot=snpAnnot1, scanAnnot=scanAnnot1)
  
  # expected output (TRUE for discordance)
  a.exp <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  a.nm <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  b.exp <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
  b.nm <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  c.exp <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
  c.nm <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
  snpID1 <- c(1,3,5,7,9)
  snpID2 <- 1:5

  # test
  a <- dp.test(genoData1, a.scanID1, snpID1,
               genoData2, a.scanID2, snpID2)
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- dp.test(genoData1, b.scanID1, snpID1,
               genoData2, b.scanID2, snpID2)
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
  c <- dp.test(genoData1, c.scanID1, snpID1,
               genoData2, c.scanID2, snpID2)
  checkIdentical(c$discordant, c.exp)
  checkIdentical(c$nonmissing, c.nm)
}



test_duplicateDiscordanceAcrossDatasets <- function() {
  # test 1: only one scan per subject in both datasets
  # first set
  # snp annotation
  snpID <- 1:10
  chrom <- rep(1L, 10)
  pos <- 101:110
  rsID <- paste("rs", pos, sep="")
  alleleA <- rep("A",10)
  alleleB <- rep("G",10)
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      alleleA=alleleA, alleleB=alleleB,
                      stringsAsFactors=FALSE)
  snpAnnot1 <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:5
  subjID <- c("a","b","c","d","e")
  sex <- rep("F",5)
  scandf <- data.frame(scanID=scanID, subjID=subjID, sex=sex,
                       stringsAsFactors=FALSE)
  scanAnnot1 <- ScanAnnotationDataFrame(scandf)
  
  geno1 <- matrix(c(rep(0,5), rep(1,5),
                   rep(1,5), rep(2,5),
                   rep(2,5), rep(0,5),
                   rep(0,10),
                   rep(0,10)), ncol=5)
  mgr <- MatrixGenotypeReader(genotype=geno1, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData1 <- GenotypeData(mgr, snpAnnot=snpAnnot1, scanAnnot=scanAnnot1)

  # second set
  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- as.integer(c(101,103,105,107,109))
  rsID <- paste("rs", pos, sep="")
  alleleA <- rep("A",5)
  alleleB <- rep("G",5)
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      alleleA=alleleA, alleleB=alleleB,
                      stringsAsFactors=FALSE)
  snpAnnot2 <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:5
  subjID <- c("c","f","b","a",NA)
  sex <- rep("F",5)
  scandf <- data.frame(scanID=scanID, subjID=subjID, sex=sex,
                       stringsAsFactors=FALSE)
  scanAnnot2 <- ScanAnnotationDataFrame(scandf)

  geno2 <- matrix(c(2,NA,1,1,0,
                   rep(0,5),
                   1,1,0,2,NA,
                   2,0,0,1,NA,
                   rep(0,5)), ncol=5)
  mgr <- MatrixGenotypeReader(genotype=geno2, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData2 <- GenotypeData(mgr, snpAnnot=snpAnnot2, scanAnnot=scanAnnot2)

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

  # expected output for minor.allele.only
  snp.exp <- data.frame(discordant=c(1,0,2,1,0), npair=c(2,1,2,2,0),
    n.disc.subj=c(1,0,2,1,0), discord.rate=c(1,0,2,1,0)/c(2,1,2,2,0))
  row.names(snp.exp) <- c("rs101","rs103","rs105","rs107","rs109")
  subj.exp <- list(c=matrix(2/2, 1, 1, dimnames=list(1,3)),
                   b=matrix(1/3, 1, 1, dimnames=list(3,2)),
                   a=matrix(1/2, 1, 1, dimnames=list(4,1)))
  # switch order to get allele freq from genoData2
  discord <- duplicateDiscordanceAcrossDatasets(genoData2, genoData1,
    rep("subjID", 2), rep("rsID", 2), minor.allele.only=TRUE)
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)
  
  # test 2: add a second scan for subject "c" to dataset 2
  scanID <- 1:5
  subjID <- c("c","f","b","a","c")
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot2 <- ScanAnnotationDataFrame(scandf)
  geno2 <- matrix(c(2,1,1,1,0,
                   rep(0,5),
                   1,1,0,2,2,
                   2,0,0,1,1,
                   2,2,2,1,0), ncol=5)
  mgr <- MatrixGenotypeReader(genotype=geno2, snpID=snpAnnot2$snpID,
    chromosome=snpAnnot2$chrom, position=snpAnnot2$pos, scanID=scanID)
  genoData2 <- GenotypeData(mgr, snpAnnot=snpAnnot2, scanAnnot=scanAnnot2)

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
    rep("subjID", 2), rep("rsID", 2), one.pair.per.subj=FALSE)
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)

  
  # check only one scan per subject
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2), one.pair.per.subj=TRUE)
  checkIdentical(discord$discordance.by.snp$discordant, discord$discordance.by.snp$n.disc.subj)
  checkEquals(3, max(discord$discordance.by.snp$npair))

  
  # select particular snps to include
  snp.include <- c("rs101","rs103","rs105")
  snp.exp <- data.frame(discordant=c(1,1,2), npair=c(4,4,4),
    n.disc.subj=c(1,1,2), discord.rate=c(1,1,2)/c(4,4,4))
  row.names(snp.exp) <- c("rs101","rs103","rs105")
  subj.exp <- list(a=matrix(1/3, 1, 1, dimnames=list(1,4)),
                   b=matrix(1/3, 1, 1, dimnames=list(2,3)),
                   c=matrix(c(2/3, 0), 1, 2, dimnames=list(3,c(1,5))))
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2), snp.include=snp.include, one.pair.per.subj=FALSE)
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)

  
  # change the order of the snps
  snp.include <- c("rs103","rs105","rs101")
  snp.exp <- data.frame(discordant=c(1,2,1), npair=c(4,4,4),
    n.disc.subj=c(1,2,1), discord.rate=c(1,2,1)/c(4,4,4))
  row.names(snp.exp) <- c("rs103","rs105","rs101")
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    rep("subjID", 2), rep("rsID", 2), snp.include=snp.include, one.pair.per.subj=FALSE)
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
    rep("subjID", 2), rep("rsID", 2), scan.exclude1=excl1, scan.exclude2=excl2,
                                                one.pair.per.subj=FALSE)
  checkIdentical(discord$discordance.by.snp, snp.exp)
  checkIdentical(discord$discordance.by.subject, subj.exp)
  
  
  # error check - snps not in dataset 2
  snp.include <- c("rs102", "rs104")
  checkException({
    discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
      rep("subjID", 2), rep("rsID", 2), snp.include=snp.include, one.pair.per.subj=FALSE)
  })

  
  # error check - datasets with no common snps (expect NULL)
  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- 111:115
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot2 <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:4
  subjID <- c("c","f","b","a")
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot2 <- ScanAnnotationDataFrame(scandf)

  geno2 <- matrix(c(rep(0,5), rep(0,5), rep(0,5), rep(0,5)), ncol=4)
  mgr <- MatrixGenotypeReader(genotype=geno2, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData2 <- GenotypeData(mgr, snpAnnot=snpAnnot2, scanAnnot=scanAnnot2)
  
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
      rep("subjID", 2), rep("rsID", 2), one.pair.per.subj=FALSE)
  checkIdentical(discord, NULL)

  
  # error check - datasets with no common subjects (expect NULL)
  # snp annotation
  snpID <- 1:10
  chrom <- rep(1L, 10)
  pos <- 101:110
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot1 <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:5
  subjID <- c("a","b","c","d",NA)
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot1 <- ScanAnnotationDataFrame(scandf)
  
  geno1 <- matrix(c(rep(0,5), rep(1,5),
                   rep(1,5), rep(2,5),
                   rep(2,5), rep(0,5),
                   rep(0,10),
                   rep(0,10)), ncol=5)
  mgr <- MatrixGenotypeReader(genotype=geno1, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData1 <- GenotypeData(mgr, snpAnnot=snpAnnot1, scanAnnot=scanAnnot1)

  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- as.integer(c(101,103,105,107,109))
  rsID <- paste("rs", pos, sep="")
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                      stringsAsFactors=FALSE)
  snpAnnot2 <- SnpAnnotationDataFrame(snpdf)

  # scan annotation
  scanID <- 1:4
  subjID <- c("g","f",NA,NA)
  scandf <- data.frame(scanID=scanID, subjID=subjID,
                       stringsAsFactors=FALSE)
  scanAnnot2 <- ScanAnnotationDataFrame(scandf)

  # netCDF
  geno2 <- matrix(c(rep(0,5), rep(0,5), rep(0,5), rep(0,5)), ncol=4)
  mgr <- MatrixGenotypeReader(genotype=geno2, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData2 <- GenotypeData(mgr, snpAnnot=snpAnnot2, scanAnnot=scanAnnot2)
  
  discord <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
      rep("subjID", 2), rep("rsID", 2), one.pair.per.subj=FALSE)
  checkIdentical(discord, NULL)
}


test_missing.fail <- function() {
  # snp annotation
  snpID <- 1:5
  chrom <- rep(1L, 5)
  pos <- 101:105
  alleleA <- rep("A",5)
  alleleB <- rep("G",5)
  snpdf <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      alleleA=alleleA, alleleB=alleleB,
                      stringsAsFactors=FALSE)
  snpAnnot <- SnpAnnotationDataFrame(snpdf)
  
  # scan annotation
  scanID <- 1:2
  subjID <- c("a","b")
  scandf <- data.frame(scanID=scanID, subjID=subjID)
  scanAnnot <- ScanAnnotationDataFrame(scandf)
  
  geno1 <- matrix(c(0,1,
                    0,1,
                    0,1,
                    0,0,
                    NA,NA), ncol=2, byrow=TRUE)
  mgr <- MatrixGenotypeReader(genotype=geno1, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData1 <- GenotypeData(mgr, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  geno2 <- matrix(c(NA,NA,
                    0,1,
                    1,0,
                    NA,NA, 
                    0,0),
                  ncol=2, byrow=TRUE)
  mgr <- MatrixGenotypeReader(genotype=geno2, snpID=snpID,
    chromosome=chrom, position=pos, scanID=scanID)
  genoData2 <- GenotypeData(mgr, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

  
  # expected output (TRUE for discordance)
  # samples "a" and "b" - default
  exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  nm <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
  a <- dp.test(genoData1, 1, 1:5,
                       genoData2, 1, 1:5)
  checkIdentical(a$discordant, exp)
  checkIdentical(a$nonmissing, nm)
  b <- dp.test(genoData1, 2, 1:5,
                      genoData2, 2, 1:5)
  checkIdentical(b$discordant, exp)
  checkIdentical(b$nonmissing, nm)

  # missing.fail[1]=TRUE
  exp <- c(FALSE, FALSE, TRUE, FALSE, TRUE)
  nm <- c(FALSE, TRUE, TRUE, FALSE, TRUE)
  a <- dp.test(genoData1, 1, 1:5,
                       genoData2, 1, 1:5,
                      missing.fail=c(TRUE,FALSE))
  checkIdentical(a$discordant, exp)
  checkIdentical(a$nonmissing, nm)
  b <- dp.test(genoData1, 2, 1:5,
                      genoData2, 2, 1:5,
                      missing.fail=c(TRUE,FALSE))
  checkIdentical(b$discordant, exp)
  checkIdentical(b$nonmissing, nm)
  
  # missing.fail[2]=TRUE
  exp <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  nm <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
  a <- dp.test(genoData1, 1, 1:5,
                       genoData2, 1, 1:5,
                      missing.fail=c(FALSE,TRUE))
  checkIdentical(a$discordant, exp)
  checkIdentical(a$nonmissing, nm)
  b <- dp.test(genoData1, 2, 1:5,
                      genoData2, 2, 1:5,
                      missing.fail=c(FALSE,TRUE))
  checkIdentical(b$discordant, exp)
  checkIdentical(b$nonmissing, nm)

  # both TRUE
  exp <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
  nm <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  a <- dp.test(genoData1, 1, 1:5,
                       genoData2, 1, 1:5,
                      missing.fail=c(TRUE,TRUE))
  checkIdentical(a$discordant, exp)
  checkIdentical(a$nonmissing, nm)
  b <- dp.test(genoData1, 2, 1:5,
                      genoData2, 2, 1:5,
                      missing.fail=c(TRUE,TRUE))
  checkIdentical(b$discordant, exp)
  checkIdentical(b$nonmissing, nm)
  
  # missing.fail[1]=TRUE and minor.allele.only=TRUE
  a.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  a.nm <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  b.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  b.nm <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
  a <- dp.test(genoData1, 1, 1:5,
                      genoData2, 1, 1:5,
                      major.genotype="G/G",
                      missing.fail=c(TRUE,FALSE))
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- dp.test(genoData1, 2, 1:5,
                      genoData2, 2, 1:5,
                      major.genotype="G/G",
                      missing.fail=c(TRUE,FALSE))
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
  
  # missing.fail[2]=TRUE and minor.allele.only=TRUE
  a.exp <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  a.nm <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  b.exp <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
  b.nm <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
  a <- dp.test(genoData1, 1, 1:5,
                      genoData2, 1, 1:5,
                      major.genotype="G/G",
                      missing.fail=c(FALSE,TRUE))
  checkIdentical(a$discordant, a.exp)
  checkIdentical(a$nonmissing, a.nm)
  b <- dp.test(genoData1, 2, 1:5,
                      genoData2, 2, 1:5,
                      major.genotype="G/G",
                      missing.fail=c(FALSE,TRUE))
  checkIdentical(b$discordant, b.exp)
  checkIdentical(b$nonmissing, b.nm)
}
