# find discordant genotype rates between duplicate scans of the same subject
# in multiple datasets

# discordantPair
# inputs:
# - two GenotypeData objects
# - two scanIDs
# - two vectors of matched snp IDs
# returns: logical vector of discordances for all snps
discordantPair <- function(genoData1, scanID1, snpID1,
                             genoData2, scanID2, snpID2) {
  # check that snp ID vectors are the same length
  stopifnot(length(snpID1) == length(snpID2))
  
  # find index of scanID1
  scanIndex1 <- which(getScanID(genoData1) == scanID1)
  if (length(scanIndex1) == 0) stop("scanID1 not found in genoData1")
  # get genotypes for this index
  geno1 <- getGenotype(genoData1, snp=c(1,-1), scan=c(scanIndex1,1))
  # discard Y chrom SNPs for females
  if (hasSex(genoData1)) {
    sex <- getSex(genoData1, index=scanIndex1)
    if (sex == "F") {
      geno1[getChromosome(genoData1, char=TRUE) == "Y"] <- NA
    }
  }
  # get matched snps
  snpIndex1 <- match(snpID1, getSnpID(genoData1))
  if (any(is.na(snpIndex1))) stop ("some SNPs not found in genoData1") 
  geno1 <- geno1[snpIndex1]
  
  # find index of scanID2
  scanIndex2 <- which(getScanID(genoData2) == scanID2)
  if (length(scanIndex2) == 0) stop("scanID2 not found in genoData2")
  # get genotypes for this index
  geno2 <- getGenotype(genoData2, snp=c(1,-1), scan=c(scanIndex2,1))
  # only need to set one set of Y genotypes to NA to avoid counting them
  # get matched snps
  snpIndex2 <- match(snpID2, getSnpID(genoData2))
  if (any(is.na(snpIndex2))) stop ("some SNPs not found in genoData2") 
  geno2 <- geno2[snpIndex2]

  # check that genotypes are the same length after selecting snps
  stopifnot(length(geno1) == length(geno2))
  
  # compare genotypes
  nonmissing <- !is.na(geno1) & !is.na(geno2)
  discordant <- nonmissing & geno1 != geno2
  return(data.frame(discordant=discordant, nonmissing=nonmissing))
}


# duplicateDiscordanceAcrossDatasets
# inputs:
# - list of GenotypeData objects
# - vector of common subject ID columns
# - vector of common snp ID columns
# - vectors of scans to exclude (optional)
# - vector of snp IDs to include (optional)
# returns: vector of number of discordances with names = common snp ID
duplicateDiscordanceAcrossDatasets <- function(genoData1, genoData2,
                                                  subjName.cols,
                                                  snpName.cols,
                                 one.pair.per.subj = TRUE,
                                               scan.exclude1=NULL,scan.exclude2=NULL,
                                                  snp.include = NULL, verbose=TRUE) {
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

  sample.annotation <- data.frame(scanID = c(scanID1, scanID2), subjID = c(subjID1, subjID2),
                                  dataset = c( rep(1, length(scanID1)), rep(2, length(scanID2))),
                                  stringsAsFactors=FALSE)
  
  dups <- intersect(sample.annotation[sample.annotation$dataset == 1, "subjID"],
                    sample.annotation[sample.annotation$dataset == 2, "subjID"])
  dups <- dups[!is.na(dups)]
  if (length(dups) == 0) {
    warning("no duplicate subjects found")
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
  
  # get snp names
  snp1 <- getSnpVariable(genoData1, snpName.cols[1])
  snp2 <- getSnpVariable(genoData2, snpName.cols[2])
  
  # find snps common to both datasets if snp.include = NULL
  if (is.null(snp.include)) {
    snp.include <- intersect(snp1, snp2)
    if (length(snp.include) == 0) {
      warning("no common snps found")
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

  nsnp <- length(snp.include)
  discord <- rep(0, nsnp)
  npair <- rep(0, nsnp)
  ndsubj <- rep(0, nsnp)
  fracList <- list(length = length(ids))
    
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
        res <- discordantPair(genoData1, scan1[i], snpID1,
                              genoData2, scan2[j], snpID2)
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
  row.names(snp.res) <- snp.include
  
  discord.res <- list()
  discord.res$discordance.by.snp <- snp.res
  discord.res$discordance.by.subject <- fracList
  return(discord.res)
}
