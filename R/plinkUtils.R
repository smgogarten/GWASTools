# ped file format is: follwing values per line,
# Family ID
# Individual ID
# Father ID
# Mother ID
# sex (1=male, 2=female, other=unknown)
# phenotype (-9=missing, 0=missing, 1=unaffected, 2=affected) 
# genotypes (both alleles e.g. A A for all SNPs ordered according to the map file order

# return genotype matrix or vector in plink format
getPlinkGenotype <- function(genoData, scan.start, scan.count,
                             scan.chromosome.filter = NULL,
                             alleleA.col=NULL, alleleB.col=NULL) {
  geno <- getGenotype(genoData, snp=c(1,-1), scan=c(scan.start,scan.count))
  if (length(dim(geno)) < 2) {
    geno <- matrix(geno, ncol=1)
  }
  # apply scan-chromosome filter
  if (!is.null(scan.chromosome.filter)) {
    filt <- scan.chromosome.filter[scan.start:(scan.start+scan.count-1),,drop=FALSE]
    chr <- getChromosome(genoData, char=TRUE)
    for (c in colnames(scan.chromosome.filter)) {
      badscans <- which(!filt[,c])
      if (length(badscans) > 0) {
        geno[chr == c, badscans] <- NA
      }
    }
  }
  
  geno[is.na(geno)] <- "0 0"
  if (is.null(alleleA.col) | is.null(alleleB.col)) {
  # maps the -1 0 1 2 genotype coding to 00 BB AB AA coding for ped files
    geno[geno %in% 0] <- "B B"
    geno[geno %in% 1] <- "A B"
    geno[geno %in% 2] <- "A A"
  } else {
    alleles <- getSnpVariable(genoData, c(alleleA.col, alleleB.col))
    names(alleles) <- c("A","B")
    # convert to character, as pmin and pmax will not work on factors
    alleles$A <- as.character(alleles$A)
    alleles$B <- as.character(alleles$B)
    aa <- paste(alleles$A, alleles$A)
    ab <- paste(pmin(alleles$A, alleles$B), pmax(alleles$A, alleles$B)) # sorted
    bb <- paste(alleles$B, alleles$B)
    for (k in 1:ncol(geno)) {
      geno[geno[,k] %in% 0, k] <- bb[geno[,k] %in% 0]
      geno[geno[,k] %in% 1, k] <- ab[geno[,k] %in% 1]
      geno[geno[,k] %in% 2, k] <- aa[geno[,k] %in% 2]
    }
  }
  if (ncol(geno) == 1) {
    geno <- as.vector(geno)
  }
  return(geno)
}


# return sample data frame in plink format
getPlinkFam <- function(genoData, family.col="family", individual.col="scanID", father.col="father", mother.col="mother", phenotype.col=NULL) {
  scan.df <- getScanVariable(genoData, c(family.col, individual.col, father.col, mother.col))
  names(scan.df) <- c("family", "individual", "father", "mother")
  scan.df$sex <- 0
  sex <- getSex(genoData)
  scan.df$sex[sex %in% "M"] <- 1
  scan.df$sex[sex %in% "F"] <- 2
  if (!is.null(phenotype.col)) {
    scan.df$phenotype <- getScanVariable(genoData, phenotype.col)
  } else {
    scan.df$phenotype <- -9
  }
  return(scan.df)
}


# return map data frame in plink format
getPlinkMap <- function(genoData, rs.col="rsID", mapdist.col=NULL) {
  map.df <- getSnpVariable(genoData, c(rs.col, "position"))
  chrom <- getChromosome(genoData, char=TRUE)
  # PLINK chromosome coding
  chrom[chrom == "X"] <- 23
  chrom[chrom == "Y"] <- 24
  chrom[chrom == "XY"] <- 25
  chrom[chrom == "M"] <- 26
  chrom[chrom == "U"] <- 0
  map.df$chromosome <- as.integer(chrom)
  if (!is.null(mapdist.col)) {
    map.df$mapdist <- getSnpVariable(genoData, mapdist.col)
  } else {
    map.df$mapdist <- 0
  }
  map.df <- map.df[,c("chromosome", rs.col, "mapdist", "position")]
  return(map.df)
}


plinkWrite <- function(genoData, pedFile="testPlink",
	family.col="family", individual.col="scanID", father.col="father", mother.col="mother", phenotype.col=NULL,
        alleleA.col=NULL, alleleB.col=NULL,
	rs.col="rsID", mapdist.col=NULL,
        scan.exclude=NULL, scan.chromosome.filter=NULL, blockSize=100,
                       verbose=TRUE){

  # name of the output_file. This will create two files: output_file.ped and output_file.map
  pedfile <- paste(pedFile,"ped",sep=".")
  mapfile <- paste(pedFile,"map",sep=".")
  
  #######################################################
  # Creating map file
  map.df <- GWASTools:::getPlinkMap(genoData, rs.col, mapdist.col)
  write.table(map.df,mapfile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  # free memory
  rm(map.df)
  
  ######################################################
  # Creating ped file
  nsample <- nscan(genoData)
  nsnp <- nsnp(genoData)
  
  scan.df <- GWASTools:::getPlinkFam(genoData, family.col, individual.col, father.col, mother.col, phenotype.col)
  scanID <- getScanID(genoData)
  
  # This loop retrives each individuals genotypes from the .nc files and writes them in the ped file. One individual per line
  nBlock <- ceiling(nsample / blockSize)
  lastB <- nsample %% blockSize
  i <- 1
  app <- FALSE
  for(j in 1:nBlock){
    bsize <- blockSize
    if(j==nBlock) bsize <- lastB
    if (verbose) message("writing block ",j," of ",nBlock,": scans ",i,"-",i+bsize-1)
    rdf <- matrix(nrow=bsize, ncol=6+nsnp)
    rdf[,1:6] <- as.matrix(scan.df[i:(i+bsize-1),])
    geno <- GWASTools:::getPlinkGenotype(genoData, scan.start=i, scan.count=bsize,
                             scan.chromosome.filter=scan.chromosome.filter,
                             alleleA.col=alleleA.col, alleleB.col=alleleB.col)
    # transpose since PLINK has scans as rows and snps as columns
    rdf[,7:(6+nsnp)] <- t(geno)
    # exclude scans
    keep <- !(scanID[i:(i+bsize-1)] %in% scan.exclude)
    rdf <- rdf[keep,]
    write.table(rdf, pedfile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=app)
    app <- TRUE
    i <- i + bsize
  }
  
  return(invisible(NULL))
}


plinkCheck <- function(genoData, pedFile, logFile="plinkCheck.txt",
	family.col="family", individual.col="scanID", father.col="father", mother.col="mother", phenotype.col=NULL,
        alleleA.col=NULL, alleleB.col=NULL,
	rs.col="rsID", map.alt=NULL,
                       check.parents=TRUE, check.sex=TRUE, 
                       scan.exclude=NULL, scan.chromosome.filter=NULL, verbose=TRUE) { 

  # return value - set to FALSE if an error is encountered
  retval <- TRUE

  # open log file
  con <- file(logFile, "w")
  
  # name of the plink files
  pedfile <- paste(pedFile,"ped",sep=".")
  mapfile <- paste(pedFile,"map",sep=".")

  # check SNPs
  if (verbose) message("Checking SNPs against map file")
  writeLines("Checking SNPs against map file...", con)
  map <- read.table(mapfile, as.is=TRUE, comment.char="")
  if (ncol(map) > 3) map <- map[,c(1,2,4)] # skip map distance
  names(map) <- c("chromosome", "rsID", "position")
  # convert chromosome to integer (PLINK coding)
  map$chromosome[map$chromosome == "X"] <- 23
  map$chromosome[map$chromosome == "Y"] <- 24
  map$chromosome[map$chromosome == "XY"] <- 25
  map$chromosome[map$chromosome == "MT"] <- 26
  snp.plink <- paste(map$chromosome, map$rsID, map$position)

  if (is.null(map.alt)) {
    map.df <- GWASTools:::getPlinkMap(genoData, rs.col=rs.col)
    map.df <- map.df[,c("chromosome", rs.col, "position")]
    names(map.df) <- c("chromosome", "rsID", "position")
  } else {
    snpID <- getSnpID(genoData)
    map.df <- map.alt[match(snpID, map.alt$snpID),
                      c("chromosome", "rsID", "position")]
  }
  snp.ncdf <- paste(map.df$chromosome, map.df$rsID, map.df$position)

  mismatch.plink <- setdiff(snp.plink, snp.ncdf)
  mismatch.ncdf <- setdiff(snp.ncdf, snp.plink)
  if (length(c(mismatch.plink, mismatch.ncdf)) > 0 ) {
    p.df <- cbind("file"=rep("ped", length(mismatch.plink)), "snp"=mismatch.plink)
    n.df <- cbind("file"=rep("netcdf", length(mismatch.ncdf)), "snp"=mismatch.ncdf)    
    mismatch <- rbind(p.df, n.df)
    writeLines("SNPs are not identical", con)
    write.table(mismatch, con, quote=FALSE, col.names=FALSE, row.names=FALSE)
    retval <- FALSE
  } else {
    writeLines("OK", con)
  }

  snp.match <- match(snp.plink, snp.ncdf)

  # sample data
  scan.df <- GWASTools:::getPlinkFam(genoData, family.col, individual.col, father.col, mother.col, phenotype.col)
  # exclude some samples
  scanID <- getScanID(genoData)
  keep <- !(scanID %in% scan.exclude)

  # read each line of PLINK file and compare with netcdf
  if (verbose) message("Checking sample data and genotypes in each line of ped file")
  writeLines("Checking sample data and genotypes in each line of ped file...", con)
  ped <- file(pedfile, "r")
  on.exit(close(ped))
  line <- 1
  # keep track of the scans we've found
  scanlist <- vector()
  while(TRUE) {
    x <- scan(ped, what="", nlines=1, quiet=TRUE)
    if (length(x) == 0) break
    if (verbose & line%%100 == 0) message("checking line ",line)

    # find matching sample in genoData
    ind <- which(scan.df[,"individual"] == x[2] & keep)
    scanlist <- c(scanlist, x[2])
    if (length(ind) == 0) {
      writeLines(paste("sample", x[2], "at line", line, "not found in NetCDF"), con)
      retval <- FALSE
      next
    } else if (length(ind) > 1) {
      writeLines(paste("sample", x[2], "at line", line, "has multiple entries in NetCDF"), con)
      retval <- FALSE
      next
    }
    
    # compare sample data
    if (!allequal(scan.df[ind,"family"], x[1])) {
      writeLines(c(paste("family mismatch for sample", x[2], "at line", line),
                   paste("Ped:", paste(x[1], collapse=" ")),
                   paste("NetCDF:", paste(scan.df[ind,"family"], collapse=" "))), con)
      retval <- FALSE
    }
    if (check.parents) {
      if (!allequal(scan.df[ind,c("father","mother")], x[3:4])) {
        writeLines(c(paste("parent mismatch for sample", x[2], "at line", line),
                     paste("Ped:", paste(x[3:4], collapse=" ")),
                     paste("NetCDF:", paste(scan.df[ind,c("father","mother")], collapse=" "))), con)
        retval <- FALSE
      }
    }
    if (check.sex) {
      if (!allequal(scan.df[ind,"sex"], x[5])) {
        writeLines(c(paste("sex mismatch for sample", x[2], "at line", line),
                     paste("Ped:", paste(x[5], collapse=" ")),
                     paste("NetCDF:", paste(scan.df[ind,"sex"], collapse=" "))), con)
        retval <- FALSE
      }
    }
    if (!is.null(phenotype.col)) {
      if (!allequal(scan.df[ind,"phenotype"], x[6])) {
        writeLines(c(paste("phenotype mismatch for sample", x[2], "at line", line),
                     paste("Ped:", x[6]),
                     paste("NetCDF:", scan.df[ind,"phenotype"])), con)
        retval <- FALSE
      }
    }

    # compare genotypes
    geno <- GWASTools:::getPlinkGenotype(genoData, scan.start=ind, scan.count=1,
                             scan.chromosome.filter=scan.chromosome.filter,
                             alleleA.col=alleleA.col, alleleB.col=alleleB.col)
    geno <- geno[snp.match]

    # sort allele by character
    a <- x[seq(7,length(x),2)]
    b <- x[seq(8,length(x),2)]
    geno.plink <- paste(pmin(a,b), pmax(a,b))
    # only check genotypes with matching SNPs
    if (length(mismatch.plink > 0)) geno.plink[snp.plink %in% mismatch.plink] <- NA
    if (length(geno.plink) != length(geno)) {
      writeLines(c("numbers of SNPs do not match",
                 paste("Ped:", length(geno.plink)),
                 paste("NetCDF:", length(geno))), con)
      retval <- FALSE
    }
    if (!allequal(geno, geno.plink)) {
      badind <- which(geno != geno.plink)
      writeLines(c(paste("genotype mismatch for sample", x[2], "at line", line, " snp", badind[1]),
                   paste("Ped:", paste(head(geno.plink[badind]), collapse=" "), "..."),
                   paste("NetCDF:", paste(head(geno[badind]), collapse=" "), "...")), con)
      retval <- FALSE
    }
    line <- line + 1
  }
  if (retval) writeLines("OK", con)

  # check that all scans were found
  writeLines("Checking that all samples were found in ped file...", con)
  missing <- setdiff(scan.df[keep,"individual"], scanlist)
  if (length(missing > 0)) {
    writeLines("samples not found in Ped:", con)
    for (i in missing) {
      writeLines(as.character(i), con)
    }
    retval <- FALSE
  } else {
    writeLines("OK", con)
  }

  close(con)
  return(retval)
}


plinkToNcdf <- function(pedFile, mapFile, nSamples,
                        ncdfFile, snpAnnotFile, scanAnnotFile,
                        ncdfXchromCode=23, ncdfXYchromCode=24, ncdfYchromCode=25,
                        ncdfMchromCode=26, ncdfUchromCode=27,
                        pedMissingCode=0, verbose=TRUE) {
  # read map file to get SNP annotation
  map <- read.table(mapFile, as.is=TRUE, comment.char="") 
  names(map)[1:4] <- c("chromosome", "rsID", "mapdist", "position")

  # are chromosomes in map file characters or integers?
  if (is.character(map$chromosome)) {
    map$chromosome[map$chromosome == "X"] <- ncdfXchromCode
    map$chromosome[map$chromosome == "XY"] <- ncdfXYchromCode
    map$chromosome[map$chromosome == "Y"] <- ncdfYchromCode
    map$chromosome[map$chromosome == "MT"] <- ncdfMchromCode
  } else {
    # assume PLINK integer coding
    chrom <- map$chromosome
    chrom[map$chromosome == 23] <- ncdfXchromCode
    chrom[map$chromosome == 24] <- ncdfYchromCode
    chrom[map$chromosome == 25] <- ncdfXYchromCode
    chrom[map$chromosome == 26] <- ncdfMchromCode
    map$chromosome <- chrom
  }
  known.chroms <- c(1:22, ncdfXchromCode, ncdfXYchromCode, ncdfYchromCode,
                   ncdfMchromCode, ncdfUchromCode)
  map$chromosome[!(map$chromosome %in% known.chroms)] <- ncdfUchromCode
  map$chromosome <- as.integer(map$chromosome)
  map$position <- as.integer(map$position)

  # save current order of SNPs in plink
  snpnames <- map$rsID
  # order by chromosome and position
  map <- map[order(map$chromosome, map$position),]
  # map new order to plink order
  snpord <- match(map$rsID, snpnames)
 
  # allele conversion
  nsnp <- nrow(map)
  if (ncol(map) >= 6) {
    names(map)[5:6] <- c("alleleA", "alleleB")
    map$alleleA[map$alleleA %in% 0] <- NA
    map$alleleB[map$alleleB %in% 0] <- NA
    bad.mono <- !is.na(map$alleleA) & !is.na(map$alleleB) & map$alleleA == map$alleleB
    if (any(bad.mono)) {
      warning("monomorphic SNPs: where alleleA=alleleB, setting alleleB to NA")
      map$alleleB[bad.mono] <- NA
    }
    alleleA <- map$alleleA
    alleleB <- map$alleleB
    findAB <- FALSE
  } else {
    alleleA <- rep(NA, nsnp)
    alleleB <- rep(NA, nsnp)
    findAB <- TRUE
  }
  
  # create integer snpID
  map$snpID <- 1:nsnp
  
  # create netCDF file
  ncdfCreate(snp.annotation=map, ncdf.filename=ncdfFile,
             n.samples=nSamples, variables="genotype")
  
  # create empty scan annotation vectors
  family <- rep(NA, nSamples)
  individ <- rep(NA, nSamples)
  father <- rep(NA, nSamples)
  mother <- rep(NA, nSamples)
  sex <- rep(NA, nSamples)
  phenotype <- rep(NA, nSamples)

  # open netCDF for writing
  nc <- open.ncdf(ncdfFile, write=TRUE)
  
  # read plink file one line at a time
  # for each line, store sample data in annotation and genotype in netCDF
  ped <- file(pedFile, "r")
  on.exit(close(ped))
  line <- 0
  while(TRUE) {
    x <- scan(ped, what="", nlines=1, quiet=TRUE)
    if (length(x) == 0) break
    line <- line + 1
    if (verbose & line%%100 == 0) message("reading line ",line)

    family[line] <- x[1]
    individ[line] <- x[2]
    father[line] <- x[3]
    mother[line] <- x[4]
    sex[line] <- x[5]
    phenotype[line] <- x[6]
  
    allele1 <- x[seq(7,length(x),2)]
    allele2 <- x[seq(8,length(x),2)]
    stopifnot(length(allele1) == nsnp)
    stopifnot(length(allele2) == nsnp)
    allele1[allele1 == pedMissingCode] <- NA
    allele2[allele2 == pedMissingCode] <- NA

    # match SNP order to annotation
    allele1 <- allele1[snpord]
    allele2 <- allele2[snpord]

    # convert to AB format
    if (findAB) {
      # if alleleA is not defined, use allele1 for A
      needAlleleA <- is.na(alleleA)
      alleleA[needAlleleA] <- allele1[needAlleleA]
      # if alleleB is not defined and alleleA != allele1, use allele1 for B
      needAlleleB <- is.na(alleleB)
      new <- !is.na(alleleA) & !is.na(allele1) & alleleA != allele1
      alleleB[needAlleleB & new] <- allele1[needAlleleB & new]
      # if alleleB is not defined and alleleA != allele2, use allele2 for B
      needAlleleB <- is.na(alleleB)
      new <- !is.na(alleleA) & !is.na(allele2) & alleleA != allele2
      alleleB[needAlleleB & new] <- allele2[needAlleleB & new]
    }
    geno1 <- rep(NA, nsnp)
    geno1[allele1 == alleleA] <- "A"
    geno1[allele1 == alleleB] <- "B"
    geno2 <- rep(NA, nsnp)
    geno2[allele2 == alleleA] <- "A"
    geno2[allele2 == alleleB] <- "B"
    geno <- paste(geno1, geno2, sep="")

    # convert to number of A alleles
    nA <- rep(NA, nsnp)
    nA[geno %in% "AA"] <- 2
    nA[geno %in% c("AB","BA")] <- 1
    nA[geno %in% c("BB")] <- 0

    # add genotype to netCDF
    put.var.ncdf(nc, "genotype", vals=nA, start=c(1,line), count=c(nsnp,1))
    
  }

  # check number of samples read
  if (line != nSamples) {
    warning("Expected ",nSamples," samples but read ",line," lines in plink file")
  }
  # define scanID
  # is individ column a unique integer?
  scanID <- 1:line
  if (allequal(individ, suppressWarnings(as.integer(individ))) &
      length(individ) == length(unique(individ))) {
    scanID <- as.integer(individ[1:line])
  }    
  # add scanID to netCDF
  put.var.ncdf(nc, "sampleID", vals=scanID, start=1, count=line)
  
  close.ncdf(nc)
  
  # save scan annotation as ScanAnnotationDataFrame
  scan.df <- data.frame(individ, family, father, mother, sex, phenotype, stringsAsFactors=FALSE)
  if (nrow(scan.df) > line) scan.df <- scan.df[1:line,]
  scan.df$scanID <- scanID
  if (allequal(scanID, scan.df$individ)) {
    scan.df$individ <- NULL
  }
  sexmf <- rep(NA, line)
  sexmf[scan.df$sex == 1] <- "M"
  sexmf[scan.df$sex == 2] <- "F"
  scan.df$sex <- sexmf
  scanAnnot <- ScanAnnotationDataFrame(scan.df)
  save(scanAnnot, file=scanAnnotFile)
  
  # save snp annotation as SnpAnnotationDataFrame
  if (findAB) {
    map$alleleA <- alleleA
    map$alleleB <- alleleB
  }
  snpAnnot <- SnpAnnotationDataFrame(map,
    XchromCode=as.integer(ncdfXchromCode), YchromCode=as.integer(ncdfYchromCode),
    XYchromCode=as.integer(ncdfXYchromCode), MchromCode=as.integer(ncdfMchromCode))
  save(snpAnnot, file=snpAnnotFile)
}
