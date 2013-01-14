# find discordant genotype rates between duplicate scans of the same subject
# in multiple datasets

.duplicatePairs <- function(genoData1, genoData2, subjName.cols,
                            scan.exclude1=NULL, scan.exclude2=NULL,
                            one.pair.per.subj=TRUE) {

  scanID1 <- getScanID(genoData1)
  subjID1 <- getScanVariable(genoData1, subjName.cols[1])
  if (!is.null(scan.exclude1)) {
    excl <- is.element(scanID1, scan.exclude1)
    scanID1 <- scanID1[!excl]
    subjID1 <- subjID1[!excl]
  }

  scanID2 <- getScanID(genoData2)
  subjID2 <- getScanVariable(genoData2, subjName.cols[2])
  if (!is.null(scan.exclude2)) {
    excl <- is.element(scanID2, scan.exclude2)
    scanID2 <- scanID2[!excl]
    subjID2 <- subjID2[!excl]
  }

  sample.annotation <- data.frame(scanID=c(scanID1, scanID2), subjID=c(subjID1, subjID2),
                                  dataset=c( rep(1, length(scanID1)), rep(2, length(scanID2))),
                                  stringsAsFactors=FALSE)
  
  dups <- intersect(sample.annotation[sample.annotation$dataset == 1, "subjID"],
                    sample.annotation[sample.annotation$dataset == 2, "subjID"])
  dups <- dups[!is.na(dups)]
  if (length(dups) == 0) {
    return(NULL)
  }
   
  ids <- list()
  for (i in 1:length(dups)) {
    ids[[i]] <- sample.annotation[is.element(sample.annotation[,"subjID"], dups[i]),
                                  c("scanID", "dataset")]
    # if one pair per subj, randomly select one sample from each dataset
    if (one.pair.per.subj) {
      for (ds in 1:2) {
        if (sum(ids[[i]]$dataset == ds) > 1) {
          ind <- which(ids[[i]]$dataset == ds)
          keep <- sample(ind, 1)
          rm <- setdiff(ind, keep)
          ids[[i]] <- ids[[i]][-rm,]
        }
      }
    }
  }
  names(ids) <- dups
  return(ids)
}

.commonSnps <- function(genoData1, genoData2, snpName.cols,
                        snp.include=NULL) {

  # get snp names
  snp1 <- getSnpVariable(genoData1, snpName.cols[1])
  snp2 <- getSnpVariable(genoData2, snpName.cols[2])
  
  # find snps common to both datasets if snp.include=NULL
  if (is.null(snp.include)) {
    snp.include <- intersect(snp1, snp2)
    if (length(snp.include) == 0) {
      return(NULL)
    }
  }
  
  # find snpIDs of common snps in each dataset
  snpID1 <- getSnpID(genoData1, match(snp.include, snp1))
  if (length(snpID1) < length(snp.include)) {
    stop("some selected snps not found in genoData1")
  }
  snpID2 <- getSnpID(genoData2, match(snp.include, snp2))
  if (length(snpID2) < length(snp.include)) {
    stop("some selected snps not found in genoData2")
  }

  snps <- data.frame(snpID1, snpID2)
  row.names(snps) <- snp.include
  return(snps)
}

.majorGenotype <- function(genoData, scanID, snpID) {
                           
  scan.excl <- setdiff(getScanID(genoData), scanID)
  allele.freq <- alleleFrequency(genoData, scan.exclude=scan.excl, verbose=FALSE)
  allele.freq <- allele.freq[as.character(snpID),"all"]

  snpIndex <- match(snpID, getSnpID(genoData))
  alleleA <- getAlleleA(genoData, index=snpIndex)
  alleleB <- getAlleleB(genoData, index=snpIndex)
  if (is.null(alleleA) | is.null(alleleB)) {
    alleleA <- rep("A", length(allele.freq))
    alleleB <- rep("B", length(allele.freq))
  }

  # find the genotype to be ignored (no minor allele) for each SNP
  major.genotype <- rep(NA, length(allele.freq))
  # A allele freq < 0.5, so A is minor allele, so BB is ignored
  Amin <- !is.na(allele.freq) & allele.freq < 0.5
  major.genotype[Amin] <- paste(alleleB[Amin], alleleB[Amin], sep="/")
  # A allele freq > 0.5, so B is minor allele, so AA is ignored
  Bmin <- !is.na(allele.freq) & allele.freq >= 0.5
  major.genotype[Bmin] <- paste(alleleA[Bmin], alleleA[Bmin], sep="/")
  return(major.genotype)
}

.selectGenotype <- function(genoData, scanID, snpID) {  
  # find index of scanID
  scanIndex <- which(getScanID(genoData) == scanID)
  if (length(scanIndex) == 0) stop("scanID not found in genoData")
  # get genotypes for this index
  geno <- getGenotype(genoData, snp=c(1,-1), scan=c(scanIndex, 1), char=TRUE)
  # discard Y chrom SNPs for females
  if (hasSex(genoData)) {
    sex <- getSex(genoData, index=scanIndex)
    if (!is.na(sex) & sex == "F") {
      geno[getChromosome(genoData, char=TRUE) == "Y"] <- NA
    }
  }
  # get matched snps
  snpIndex <- match(snpID, getSnpID(genoData))
  if (any(is.na(snpIndex))) stop ("some SNPs not found in genoData") 
  geno <- geno[snpIndex]
  return(geno)
}

# discordantPair
# inputs: two genotype vectors
# returns: logical vector of discordances for all snps
# if major.genotype is not NULL, "nonmissing" only counts SNPs
#  where the pair included the minor allele
.discordantPair <- function(geno1, geno2,
                            major.genotype=NULL,
                            missing.fail=c(FALSE, FALSE)) {
  stopifnot(length(geno1) == length(geno2))
  
  # compare genotypes
  nonmissing <- !is.na(geno1) & !is.na(geno2)
  if (!is.null(major.genotype)) {
    nonmissing <- nonmissing & !is.na(major.genotype) & (geno1 != major.genotype | geno2 != major.genotype)
  }
  discordant <- nonmissing & geno1 != geno2
  if (missing.fail[1]) {
    fail <- is.na(geno1) & !is.na(geno2)
    if (!is.null(major.genotype)) {
      fail <- fail & !is.na(major.genotype) & (geno2 != major.genotype)
    }
    discordant <- discordant | fail
    nonmissing <- nonmissing | fail
  }
  if (missing.fail[2]) {
    fail <- !is.na(geno1) & is.na(geno2)
    if (!is.null(major.genotype)) {
      fail <- fail & !is.na(major.genotype) & (geno1 != major.genotype)
    }
    discordant <- discordant | fail
    nonmissing <- nonmissing | fail
  }
  return(data.frame(discordant=discordant, nonmissing=nonmissing))
}

# duplicateDiscordanceAcrossDatasets
# inputs:
# - list of GenotypeData objects
# - vector of common subject ID columns
# - vector of common snp ID columns
# - vectors of scans to exclude (optional)
# - vector of snp IDs to include (optional)
duplicateDiscordanceAcrossDatasets <- function(genoData1, genoData2,
                                               subjName.cols, snpName.cols,
                                               one.pair.per.subj=TRUE,
                                               minor.allele.only=FALSE,
                                               missing.fail=c(FALSE,FALSE),
                                               scan.exclude1=NULL,scan.exclude2=NULL,
                                               snp.include=NULL, verbose=TRUE) {
  # check that both genoData objects have subjName, snpName
  stopifnot(hasScanVariable(genoData1, subjName.cols[1]))
  stopifnot(hasSnpVariable(genoData1, snpName.cols[1]))
  stopifnot(hasScanVariable(genoData2, subjName.cols[2]))
  stopifnot(hasSnpVariable(genoData2, snpName.cols[2]))

  if ("Y" %in% c(getChromosome(genoData1, char=TRUE))) {
    if (!hasSex(genoData1)) {
      stop("sex is required for checking Y chromosome discordance")
    }
  }
  
  # find duplicate scans
  ids <- .duplicatePairs(genoData1, genoData2, subjName.cols,
                         scan.exclude1, scan.exclude2,
                         one.pair.per.subj)
  if (is.null(ids)) {    
    warning("no duplicate IDs found; check subjName.cols")
    return(NULL)
  }

  # find common snps
  snps <- .commonSnps(genoData1, genoData2, snpName.cols,
                      snp.include)
  if (is.null(snps)) {
    warning("no common snps found; check snpName.cols")
    return(NULL)
  }
  
  if (minor.allele.only) {
    # calculate allele frequency of dataset with fewer snps, common samples only
    if (verbose) message("Calculating allele freqency in genoData1")

    scan.freq <- unlist(lapply(ids, function(x) {x$scanID[x$dataset == 1][1]}),
                        use.names=FALSE)
    major.genotype <- .majorGenotype(genoData1, scan.freq, snps$snpID1)
  } else {
    major.genotype <- NULL
  }
   
  nsnp <- nrow(snps)
  discord <- rep(0, nsnp)
  npair <- rep(0, nsnp)
  ndsubj <- rep(0, nsnp)
  fracList <- list(length=length(ids))
   
  # for each duplicate, calculate pair discordance
  # add to total number of discordances for each snp
  for (k in 1:(length(ids))) {
    idk <- ids[[k]] # all scanIDs for the kth dup

    n <- nrow(idk)  # total number of scans
    n1 <- sum(idk$dataset == 1) # number of scans in dataset1
    n2 <- sum(idk$dataset == 2) # number of scans in dataset2
    scan1 <- idk$scanID[idk$dataset == 1]
    scan2 <- idk$scanID[idk$dataset == 2]
    
    if (verbose)  
      message("subject ",k, " out of ",length(ids),", ",n," replications")

    frac <- matrix(NA, n1, n2, dimnames=list(scan1, scan2))
    nds <- rep(0, nsnp)
    for (i in 1:n1) {
      for (j in 1:n2) {
        # get matching genotypes
        geno1 <- .selectGenotype(genoData1, scan1[i], snps$snpID1)
        geno2 <- .selectGenotype(genoData2, scan2[j], snps$snpID2)
        res <- .discordantPair(geno1, geno2,
                               major.genotype, missing.fail)
        discord[res$discordant] <- discord[res$discordant] + 1
        npair[res$nonmissing] <- npair[res$nonmissing] + 1
        nds[res$discordant] <- nds[res$discordant] + 1
        frac[i,j] <- sum(res$discord) / sum(res$nonmissing)
      }
    }
    # discordance by snp
    nds[nds > 1] <- 1
    ndsubj <- ndsubj + nds
    
    # discordance by subject
    fracList[[k]] <- frac
  }
  names(fracList) <- names(ids)
  
  #n.disc.subj = n.subj.with.at.least.one.discordance
  snp.res <- data.frame(discordant=discord, npair=npair, n.disc.subj=ndsubj, discord.rate=discord/npair)
  row.names(snp.res) <- row.names(snps)
  
  discord.res <- list()
  discord.res$discordance.by.snp <- snp.res
  discord.res$discordance.by.subject <- fracList
  return(discord.res)
}

##########
# functions for minor allele sensitivity and specificity

.genoClass <- function(geno, major.genotype) {
  a <- substr(geno, 1, 1)
  b <- substr(geno, 3, 3)
  major <- substr(major.genotype, 1, 1)
    
  class <- rep("miss", length(geno))
  class[!is.na(geno) & geno == major.genotype] <- "major"
  class[!is.na(geno) & a == major & b != major] <- "het"
  class[!is.na(geno) & a != major & b != major] <- "minor"
  return(class)
}

## test <- matrix(c("2TP", "1TP+1FP", "2FP",
##                  "1TP+1FN", "1TN+1TP", "1TN+1FP",
##                  "2FN", "1FN+1TN", "2TN",
##                  "2FN", "1FN+*", "2*"),
##                ncol=3, nrow=4, byrow=TRUE,
##                dimnames=list(c("II", "AI", "AA", "--"),
##                  c("II", "AI", "AA")))
## test
## rows: geno2, cols: geno1
##    II        AI        AA       
## II "2TP"     "1TP+1FP" "2FP"    
## AI "1TP+1FN" "1TN+1TP" "1TN+1FP"
## AA "2FN"     "1FN+1TN" "2TN"    
## -- "2FN"     "1FN+*"   "2*"   
## * = exclude from the counts
## alternatively, could treat "--" like "AA"
## or could ignore "--"
## "II"=minor, "AI"=het, "AA"=major, "--"=miss

.truePos <- function(geno1, geno2) {
  2*(geno1 == "minor" & geno2 == "minor") + 
   (geno1 == "minor" & geno2 == "het") +
   (geno1 == "het" & geno2 == "minor") +
   (geno1 == "het" & geno2 == "het")
}

.trueNeg <- function(geno1, geno2) {
  (geno1 == "het" & geno2 == "het") +
   (geno1 == "het" & geno2 == "major") +
   (geno1 == "major" & geno2 == "het") +
   2*(geno1 == "major" & geno2 == "major")
}

.falsePos <- function(geno1, geno2) {
  (geno1 == "het" & geno2 == "minor") +
   2*(geno1 == "major" & geno2 == "minor") +
   (geno1 == "major" & geno2 == "het")
}

.falseNeg <- function(geno1, geno2) {
  (geno1 == "minor" & geno2 == "het") +
   2*(geno1 == "minor" & geno2 == "major") +
   2*(geno1 == "minor" & geno2 == "miss") +
   (geno1 == "het" & geno2 == "major") +
   (geno1 == "het" & geno2 == "miss")
}


minorAlleleDetectionAccuracy <- function(genoData1, genoData2,
                                         subjName.cols, snpName.cols,
                                         scan.exclude1=NULL,scan.exclude2=NULL,
                                         snp.include=NULL, verbose=TRUE) {
  # check that both genoData objects have subjName, snpName
  stopifnot(hasScanVariable(genoData1, subjName.cols[1]))
  stopifnot(hasSnpVariable(genoData1, snpName.cols[1]))
  stopifnot(hasScanVariable(genoData2, subjName.cols[2]))
  stopifnot(hasSnpVariable(genoData2, snpName.cols[2]))

  if ("Y" %in% c(getChromosome(genoData1, char=TRUE))) {
    if (!hasSex(genoData1)) {
      stop("sex is required for checking Y chromosome")
    }
  }
  
  # find duplicate scans
  ids <- .duplicatePairs(genoData1, genoData2, subjName.cols,
                         scan.exclude1, scan.exclude2,
                         one.pair.per.subj=TRUE)
  if (is.null(ids)) {    
    warning("no duplicate IDs found; check subjName.cols")
    return(NULL)
  }

  # find common snps
  snps <- .commonSnps(genoData1, genoData2, snpName.cols,
                      snp.include)
  if (is.null(snps)) {
    warning("no common snps found; check snpName.cols")
    return(NULL)
  }
  
    # calculate allele frequency of dataset with fewer snps, common samples only
  if (verbose) message("Calculating allele freqency in genoData1")
  scan.freq <- unlist(lapply(ids, function(x) {x$scanID[x$dataset == 1][1]}),
                      use.names=FALSE)
  major.genotype <- .majorGenotype(genoData1, scan.freq, snps$snpID1)
   
  nsnp <- nrow(snps)
  truePos <- rep(0, nsnp)
  trueNeg <- rep(0, nsnp)
  falsePos <- rep(0, nsnp)
  falseNeg <- rep(0, nsnp)
  npair <- rep(0, nsnp)
   
  for (k in 1:(length(ids))) {
    idk <- ids[[k]] # all scanIDs for the kth dup
    scan1 <- idk$scanID[idk$dataset == 1]
    scan2 <- idk$scanID[idk$dataset == 2]
    
    if (verbose)  
      message("subject ",k, " out of ",length(ids))

    geno1 <- .selectGenotype(genoData1, scan1, snps$snpID1)
    geno2 <- .selectGenotype(genoData2, scan2, snps$snpID2)
    class1 <- .genoClass(geno1, major.genotype)
    class2 <- .genoClass(geno2, major.genotype)
    truePos <- truePos + .truePos(class1, class2)
    trueNeg <- trueNeg + .trueNeg(class1, class2)
    falsePos <- falsePos + .falsePos(class1, class2)
    falseNeg <- falseNeg + .falseNeg(class1, class2)
    npair <- npair + !is.na(geno1)
  }

  sensitivity <- truePos / (truePos + falseNeg)
  specificity <- trueNeg / (trueNeg + falsePos)
  positivePredictiveValue <- truePos / (truePos + falsePos)
  negativePredictiveValue <- trueNeg / (trueNeg + falseNeg)
  res <- data.frame(npair, sensitivity, specificity,
                    positivePredictiveValue, negativePredictiveValue)
  row.names(res) <- row.names(snps)
  return(res)
}
