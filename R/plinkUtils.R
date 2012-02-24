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
    geno[geno == 0] <- "B B"
    geno[geno == 1] <- "A B"
    geno[geno == 2] <- "A A"
  } else {
    alleles <- getSnpVariable(genoData, c(alleleA.col, alleleB.col))
    names(alleles) <- c("A","B")
    aa <- paste(alleles$A, alleles$A)
    ab <- paste(pmin(alleles$A, alleles$B), pmax(alleles$A, alleles$B)) # sorted
    bb <- paste(alleles$B, alleles$B)
    for (k in 1:ncol(geno)) {
      geno[geno[,k] == 0, k] <- bb[geno[,k] == 0]
      geno[geno[,k] == 1, k] <- ab[geno[,k] == 1]
      geno[geno[,k] == 2, k] <- aa[geno[,k] == 2]
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
	rs.col="rsID",
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
  map <- read.table(mapfile, as.is=TRUE, comment="")
  if (ncol(map) > 3) map <- map[,c(1,2,4)] # skip map distance
  names(map) <- c("chromosome", "rsID", "position")
  snp.plink <- paste(map$chromosome, map$rsID, map$position)
  
  map.df <- GWASTools:::getPlinkMap(genoData, rs.col=rs.col)
  map.df <- map.df[,c("chromosome", "rsID", "position")]
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
      if (scan.df[ind,"phenotype"] != x[6]) {
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
