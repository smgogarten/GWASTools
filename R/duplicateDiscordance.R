###### calculates genotype discordance for duplicate samples
# returns discordance by snp, discordance by subject, and correlation

duplicateDiscordance <- function(genoData, # object of type GenotypeData
                                 subjName.col,
                                 one.pair.per.subj = TRUE,
                                 corr.by.snp = FALSE,
                                 minor.allele.only = FALSE,
                                 allele.freq = NULL,
                                 scan.exclude = NULL,
                                 snp.exclude = NULL,
                                 snp.block.size=5000, # for by-snp correlation
                                 verbose = TRUE) {
  
  # check that column with subject IDs is included in genoData
  if (!hasScanVariable(genoData, subjName.col))
    stop(paste(subjName.col, "not found in genoData"))

  # check that sex is included for checking Y chromosome
  chrom <- getChromosome(genoData, char=TRUE)
  if ("Y" %in% chrom) {
    if (!hasSex(genoData)) {
      stop("sex is required for checking Y chromosome discordance")
    } else {
      sex <- getSex(genoData)
    }
  }

  # check that allele.freq is given if needed
  if (minor.allele.only) {
    if (is.null(allele.freq)) {
      stop("vector of A allele frequency required to determine minor allele")
    } else {
      stopifnot(is.numeric(allele.freq) & length(allele.freq) == nsnp(genoData))
    }
  }
  
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
  dups <- unique(na.omit(sample.annotation[sample.annotation$duplicated, "subjID"]))

  ids <- list()
  for (i in 1:length(dups)) {
    ids[[i]] <- sample.annotation[is.element(sample.annotation[,"subjID"], dups[i]), "scanID" ]
    # if one pair per subj, randomly select two samples
    if (one.pair.per.subj) {
      ids[[i]] <- sort(sample(ids[[i]], 2))
    }
  }
  names(ids) <- dups

  # get the indices of snps to be included
  snpID <- getSnpID(genoData)
  index <- which(!is.element(snpID, snp.exclude))

  if (minor.allele.only) {
    # find the genotype to be ignored (no minor allele) for each SNP
    major.genotype <- rep(NA,length(index))
    # A allele freq < 0.5, so A is minor allele, so BB=0 is ignored
    major.genotype[allele.freq[index] < 0.5] <- 0
    # A allele freq > 0.5, so B is minor allele, so AA=2 is ignored
    major.genotype[allele.freq[index] >= 0.5] <- 2
  }

  nsnp <- length(index)
  discord <- rep(0,nsnp)
  npair <- rep(0,nsnp)
  ndsubj <- rep(0,nsnp)
  fracList <- list(length=length(ids))
  corrList <- list(length=length(ids))

  # for each set of duplicates (which may have >3 members)
  for (k in 1:(length(ids))) {
    # get the indices of samples in the dup set
    n <- length(ids[[k]])
    idk <- which(is.element(scanID, ids[[k]]))

    if (verbose)  
      message("subject ",k, " out of ",length(ids),", ",n," replications")

    # get the genotypic data and store in dat
    dat <- matrix(NA, length(snpID), n)
    for(m in 1:n) dat[,m] <- getGenotype(genoData, snp=c(1,-1), scan=c(idk[m],1))

    # if this is a female, exclude Y chrom
    if ("Y" %in% chrom & hasSex(genoData)) {
      if (sex[idk[1]] == "F") {
        dat[chrom == "Y",] <- NA
      }
    }

    # remove snps to be excluded
    dat <- dat[index,]

    frac <- matrix(NA, n, n) # fraction of calls identical
    corr <- matrix(NA, n, n) # correlation

    nds <- rep(0,nsnp)
    for (i in 1:(n-1)) {
      nai <- !is.na(dat[,i])
      for (j in (i+1):n) {
        # naij is where both samples are non-missing
	naij <- nai & !is.na(dat[,j])
        if (minor.allele.only) {
          # naij is also where at least one sample has the minor allele
          naij <- naij & !is.na(major.genotype) & (dat[,i] != major.genotype | dat[,j] != major.genotype)
        }
	npair = npair + (naij)
        # discij is where both samples are naij AND discordant
        discij <- naij & dat[,i] != dat[,j]
	discord = discord + discij
	nds = nds | discij
	frac[i,j] <- sum(dat[,i][naij] == dat[,j][naij]) / sum(naij)
	frac[j,i] <- frac[i,j]
        # check for no variation in one of the vectors - correlation not defined in this case
        if (length(unique(dat[,i][naij])) == 1 | length(unique(dat[,j][naij])) == 1) {
          corr[i,j] <- NA
        } else {
          corr[i,j] <- cor(dat[,i][naij], dat[,j][naij])
        }
        corr[j,i] <- corr[i,j]
      }
    }
    # discordance by snp
    nds[nds > 1] <- 1
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
  
  # correlation by SNP
  if (corr.by.snp) {

    message("calculating dosage correlation by SNP, in blocks of ",
            prettyNum(snp.block.size, big.mark=","), " SNPs")
    
    # make 2 vectors with one sample from each subject (returns index)
    scan1 <- rep(NA,length(ids))
    scan2 <- rep(NA,length(ids))
    for (k in 1:(length(ids))) {
      idk <- which(is.element(scanID, ids[[k]]))
      scan1[k] <- idk[1]
      scan2[k] <- idk[2]
    }

    # make 2 matrices with one sample from each subject
    # each has all (selected) snp genotypes
    snpID.get <- snpID[index]
    geno1.matx <- .dosCorSelectGenotype(genoData, scanIDs=scanID[scan1], snpIDs=snpID.get)
    geno2.matx <- .dosCorSelectGenotype(genoData, scanIDs=scanID[scan2], snpIDs=snpID.get)

    # check snp order
    if(!allequal(rownames(geno1.matx), rownames(geno2.matx))) {
      stop("genotype matrices differ in SNP dimension")}

    # check sample order
    subj1 <- sample.annotation$subjID[match(colnames(geno1.matx), sample.annotation$scanID)]
    subj2 <- sample.annotation$subjID[match(colnames(geno2.matx), sample.annotation$scanID)]       
    if(!allequal(subj1, subj2)) stop("genotype matrices differ in sample dimension")

    # if ignoring major homozygous genotypes, set to missing where *both* samps are maj hom. 
    if (minor.allele.only) {    
        a.maj <- ifelse(major.genotype %in% 2, TRUE, FALSE)

        stopifnot(nsnp == nrow(geno1.matx))
        nscan <- length(ids)
        
        # index the matrices
        # make matrix of logical for a.maj at every genotype
        a.maj.matx <- matrix(rep(a.maj, nscan), nrow=nsnp)

        # make matrix where TRUE when both samps are major homozygous for AA
        bothAA <- a.maj.matx & geno1.matx %in% 2 & geno2.matx %in% 2

        # make matrix where TRUE when both samps are major homozygous for BB
        bothBB <- !a.maj.matx & geno1.matx %in% 0 & geno2.matx %in% 0

        # matrix of logical where TRUE when genotypes at that index should be NA
        bothHom <- bothAA | bothBB

        geno1.matx[bothHom] <- NA
        geno2.matx[bothHom] <- NA        

        ### this syntax sets to missing any major homozygous genotype:
        ## geno1.matx[a.maj, ][geno1.matx[a.maj, ] %in% 2] <- NA
        ## geno1.matx[!a.maj, ][geno1.matx[!a.maj, ] %in% 0] <- NA
        ## geno2.matx[a.maj, ][geno2.matx[a.maj, ] %in% 2] <- NA
        ## geno2.matx[!a.maj, ][geno2.matx[!a.maj, ] %in% 0] <- NA
      }

        # look through blocks of snps
        corr.snp <- rep(NA,nsnp)
        last.row <- 0
        nblocks <- ceiling(nsnp / snp.block.size)
        for (i in 1:nblocks) {

          if(verbose){message("Block ",i," of ", nblocks)}

          idx <- (1:snp.block.size) + (i - 1) * snp.block.size

          # account for where there may be less snps in the block than the block size,
          # i.e., for final block
          if(i %in% max(nblocks)) {idx <- (last.row+1):nsnp}

          # suppressWarnings will avoid the warning messages where variance SD is 0
          r.block <- suppressWarnings(diag(cor(t(geno1.matx[idx,,drop=FALSE]),
                                             t(geno2.matx[idx,,drop=FALSE]),
                                             use="pairwise.complete.obs")))

          corr.snp[idx] <- r.block
          last.row <- max(idx)

      }
    
    snp.res$correlation <- corr.snp
  }

  discord.res <- list()
  discord.res$discordance.by.snp <- snp.res
  discord.res$discordance.by.subject <- fracList
  discord.res$correlation.by.subject <- corrList
  return(discord.res)
}

