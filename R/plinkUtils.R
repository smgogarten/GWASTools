# ped file format is: follwing values per line,
# Family ID
# Individual ID
# Father ID
# Mother ID
# sex (1=male, 2=female, other=unknown)
# phenotype (-9=missing, 0=missing, 1=unaffected, 2=affected) 
# genotypes (both alleles e.g. A A for all SNPs ordered according to the map file order

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
    ab <- paste(alleles$A, alleles$B)
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

plinkWrite <- function(genoData, pedFile="testPlink", scan.chromosome.filter=NULL,
	family.col="family", individual.col="scanID", father.col="father", mother.col="mother", phenotype.col=NULL,
        alleleA.col=NULL, alleleB.col=NULL,
	rs.col="rsID", mapdist.col=NULL, scan.exclude=NULL, blockSize=100){

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



plinkCheck <- function(genoData, pedFile, scan.chromosome.filter=NULL,
	family.col="family", individual.col="scanID", father.col="father", mother.col="mother", phenotype.col=NULL,
        alleleA.col=NULL, alleleB.col=NULL,
	rs.col="rsID") { 
  
  # name of the plink files
  pedfile <- paste(pedFile,"ped",sep=".")
  mapfile <- paste(pedFile,"map",sep=".")

  # check SNPs
  map <- read.table(mapfile, as.is=TRUE, comment="")
  if (ncol(map) > 3) map <- map[,c(1,2,4)] # skip map distance
  names(map) <- c("chromosome", "rsID", "position")
  snp.plink <- paste(map$chromosome, map$rsID, map$position)
  
  map.df <- GWASTools:::getPlinkMap(genoData, rs.col=rs.col)
  map.df <- map.df[,c("chromosome", "rsID", "position")]
  snp.ncdf <- paste(map.df$chromosome, map.df$rsID, map.df$position)

  mismatch.plink <- setdiff(snp.plink, snp.ncdf)
  names(mismatch.plink) = rep("ped", length(mismatch.plink))
  mismatch.ncdf <- setdiff(snp.ncdf, snp.plink)
  names(mismatch.ncdf) = rep("netcdf", length(mismatch.ncdf))
  mismatch <- c(mismatch.plink, mismatch.ncdf)
  if (length(mismatch) > 0 ) {
    write.table(mismatch, file="snp_mismatch.log",
                quote=FALSE, col.names=FALSE)
    message("SNPs are not identical - see 'snp_mismatch.log'")
    return(FALSE)
  }

  snp.match <- match(snp.plink, snp.ncdf)

  # sample data
  scan.df <- GWASTools:::getPlinkFam(genoData, family.col, individual.col, father.col, mother.col, phenotype.col)

  # read each line of PLINK file and compare with netcdf
  ped <- file(pedfile, "r")
  on.exit(close(ped))
  line <- 1
  # keep track of the scans we've found
  scanlist <- vector()
  while(TRUE) {
    x <- scan(ped, what="", nlines=1, quiet=TRUE)
    if (length(x) == 0) break

    # find matching sample in genoData
    ind <- which(scan.df[,"individual"] == x[2])
    scanlist <- c(scanlist, ind)
    if (length(ind) == 0) {
      message("sample ", x[2], " at line ", line, " not found in NetCDF")
      return(FALSE)
    }
    
    # compare sample data
    if (!allequal(scan.df[ind,1:5], x[1:5])) {
      message("sample data mismatch for sample ", x[2], " at line ", line)
      message("Ped: ", paste(x[1:5], collapse=" "))
      message("NetCDF: ", paste(scan.df[ind,1:5], collapse=" "))
      return(FALSE)
    }
    if (!is.null(phenotype.col)) {
      if (scan.df[ind,"phenotype"] != x[6]) {
        message("phenotype mismatch for sample ", x[2], " at line ", line)
        message("Ped: ", x[6])
        message("NetCDF: ", scan.df[ind,"phenotype"])
        return(FALSE)
      }
    }

    # compare genotypes
    geno <- GWASTools:::getPlinkGenotype(genoData, scan.start=ind, scan.count=1,
                             scan.chromosome.filter=scan.chromosome.filter,
                             alleleA.col=alleleA.col, alleleB.col=alleleB.col)
    geno <- geno[snp.match]

    geno.plink <- paste(x[seq(7,length(x),2)], x[seq(8,length(x),2)])
    if (length(geno.plink) != length(geno)) {
      message("numbers of SNPs do not match")
      message("Ped: ", length(geno.plink))
      message("NetCDF: ", length(geno))
      return(FALSE)
    }
    if (!allequal(geno, geno.plink)) {
      message("genotype mismatch for sample ", x[2], " at line ", line)
      message("Ped: ", paste(head(geno.plink), collapse=" "))
      message("NetCDF: ", paste(head(geno), collapse=" "))
      return(FALSE)
    }
    line <- line + 1
  }

  # check that all scans were found
  if (length(scanlist) < length(scan.df[,"individual"])) {
    message("samples not found in Ped:")
    for (i in scan.df[-1*scanlist,"individual"]) {
      message(i)
    }
    return(FALSE)
  }
  
  return(TRUE)
}
