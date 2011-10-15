###### calculates genotype discordance for duplicate samples
# returns discordance by snp, discordance by subject, and correlation

duplicateDiscordance <- function(genoData, # object of type GenotypeData
                                  subjName.col, 
                                  scan.exclude = NULL,
                                  snp.exclude = NULL,
                                  verbose = TRUE) {
  
  # check that column with subject IDs is included in genoData
  if (!hasScanVariable(genoData, subjName.col))
    stop(paste(subjName.col, "not found in genoData"))

  # get scan and subject IDs
  scanID <- getScanID(genoData)
  subjID <- getScanVariable(genoData, subjName.col)
  sample.annotation <- data.frame("scanID"=scanID, "subjID"=subjID)
  
  if (!is.null(scan.exclude)){
    excl <- is.element(scanID, scan.exclude)
    sample.annotation <- sample.annotation[!excl,]
  }

  # find duplicate scans
  sample.annotation$duplicated = duplicated(sample.annotation[, "subjID"])
  dups <- unique(sample.annotation[sample.annotation$duplicated, "subjID"])

  ids <- list()
  for (i in 1:length(dups)) {
    ids[[i]] <- sample.annotation[ is.element(sample.annotation[,"subjID"], dups[i]), "scanID" ]
  }
  names(ids) <- dups

  # get the indices of snps to be included
  snpID <- getSnpID(genoData)
  index <- which(!is.element(snpID, snp.exclude))

  nsnp <- length(index)
  discord <- rep(0,nsnp)
  npair <- rep(0,nsnp)
  ndsubj <- rep(0,nsnp)
  fracList <- list(length=length(ids))
  corrList <- list(length=length(ids))

  # for each set of duplicates (which may have >3 members)
  for(k in 1:(length(ids))) 
  {
    if (verbose)  
      message("subject ",k, " out of ",length(ids),", ")
    
    # get the indices of samples in the dup set
    idk <-  ids[[k]]
    n <- length(idk)
    idk <- which(is.element(scanID,idk))

    if (verbose)  
      message(n," replications\n")

    # get the genotypic data and store in dat
    dat <- matrix(NA, length(snpID), n)
    for(m in 1:n) dat[,m] <- getGenotype(genoData, snp=c(1,-1), scan=c(idk[m],1))

    # remove snps to be excluded
    dat <- dat[index,]

    frac <- matrix(NA, n, n) # fraction of calls identical
    corr <- matrix(NA, n, n) # correlation

    nds <- rep(0,nsnp)
    for(i in 1:(n-1)){
      nai <- !is.na(dat[,i])
      for(j in (i+1):n){
	naij <- nai & !is.na(dat[,j])
	npair = npair  + (naij)
	discord = discord + (naij & dat[,i]!=dat[,j])
	nds = nds | (naij & dat[,i]!=dat[,j])
	frac[i,j] <- sum(dat[,i][naij]==dat[,j][naij])/sum(naij)
	frac[j,i] <- frac[i,j]
        corr[i,j] <- cor(dat[,i][naij],dat[,j][naij])
        corr[j,i] <- corr[i,j]
      }
    }
    # discordance by snp
    nds[nds>1] <- 1
    ndsubj <- ndsubj + nds
    
    # discordance by subject
    row.names(frac) <- as.character(ids[[k]])
    colnames(frac) <- as.character(ids[[k]])
    diag(frac) <- 1
    frac <- 1-frac # convert from concordance to discordance
    fracList[[k]] <- frac

    # correlation by subject
    row.names(corr) <- as.character(ids[[k]])
    colnames(corr) <- as.character(ids[[k]])
    diag(corr) <- 1
    corrList[[k]] <- corr
    
  }
  names(fracList) <- names(ids)
  names(corrList) <- names(ids)

  #n.disc.subj = n.subj.with.at.least.one.discordance
  snp.res <- data.frame(snpID=snpID[index], discordant=discord, npair=npair, n.disc.subj=ndsubj, discord.rate=discord/npair)
  
  discord.res <- list()
  discord.res$discordance.by.snp <- snp.res
  discord.res$discordance.by.subject <- fracList
  discord.res$correlation.by.subject <- corrList
  return(discord.res)
}

