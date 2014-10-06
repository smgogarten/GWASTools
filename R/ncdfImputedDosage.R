ncdfImputedDosage <- function(input.files, ncdf.filename, chromosome,
                              input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                              input.dosage=FALSE, block.size=5000,
                              snp.annot.filename="dosage.snp.RData",
                              scan.annot.filename="dosage.scan.RData",
                              verbose=TRUE) {
  
  .Deprecated("imputedDosageFile")
  
  # arguments: input type (impute2, beagle, mach), probs or dosages
  input.type <- match.arg(input.type)
  
  # determine number of SNPs and samples
  if (verbose) message("Determining number of SNPs and samples...")
  Cnt <- count.fields(input.files[1])
  line.cnt <- length(Cnt)
  col.cnt <- max(Cnt)
  if (input.type == "IMPUTE2") {
    nsnp <- line.cnt
    nsamp <- (col.cnt-5)/3
  }
  if (input.type == "BEAGLE") {
    nsnp <- line.cnt - 1
    if (input.dosage) nsamp <- col.cnt-3 else nsamp <- (col.cnt-3)/3
  }
  if (input.type == "MaCH") {
    if (input.dosage) nsnp <- col.cnt-2 else nsnp <- (col.cnt-2)/2
    nsamp <- line.cnt
  }
  snpID <- 1:nsnp
  scanID <- 1:nsamp
  
  # create NetCDF
  # define dimensions
  snpdim <- dim.def.ncdf("snp", "count", snpID)
  sampledim <- dim.def.ncdf("sample", "count", 1:nsamp, unlim=TRUE)
  chardim <- dim.def.ncdf("nchar", "", 1)

  # define variables
  varID <- var.def.ncdf("sampleID", "id", dim=sampledim, missval=0, prec="integer")
  varpos <- var.def.ncdf("position", "bases", dim=snpdim, missval=-1, prec="integer")
  varchr <- var.def.ncdf("chromosome", "id", dim=snpdim, missval=-1, prec="integer")
  varA <- var.def.ncdf("alleleA", "allele", dim=list(chardim,snpdim), missval="0", prec="char")
  varB <- var.def.ncdf("alleleB", "allele", dim=list(chardim,snpdim), missval="0", prec="char")
  vargeno <- var.def.ncdf("genotype", "A_allele_dosage", dim=list(snpdim,sampledim), missval=-1, prec="single")

  # create the NetCDF file
  if (verbose) message("Creating NetCDF file with ", nsnp, " SNPs and ", nsamp, " samples...")
  nc <- create.ncdf(ncdf.filename, list(varID, varpos, varchr, varA, varB, vargeno))
        
  # read input file(s)
  if (input.type == "IMPUTE2") {
    if (input.dosage) stop("input.dosage=TRUE not valid for input.type=IMPUTE2")
    # .samples file
    if (verbose) message ("Reading sample file...")
    samp.header <- scan(input.files[2], what=character(), nlines=1, quiet=TRUE)
    samples <- read.table(input.files[2], as.is=TRUE, header=FALSE, skip=2)
    names(samples) <- samp.header
    if (nrow(samples) != nsamp) stop("Sample number mismatch: ", nsamp, " in genotype file; ", nrow(samples), "in samples file")

    # empty matrix for SNP data
    snps <- matrix("", nrow=nsnp, ncol=5)
    
    # .gens file
    if (verbose) message ("Reading genotype file...")
    opfile <- file(input.files[1], "r")
    cnt <- 1
    while (length(ss <- scan(opfile, what=character(), quiet=TRUE, nlines=block.size)) > 0) {
      if (verbose) message("Block ", ceiling(cnt/block.size), " of ", ceiling(nsnp/block.size))
      dat <- matrix(ss, ncol=col.cnt, byrow=TRUE)
      snps[cnt:(cnt+nrow(dat)-1),] <- dat[,1:5]
      dat <- dat[,6:ncol(dat),drop=FALSE]
      mode(dat) <- "numeric"
      dosage <- .probToDosage(dat)
      if (ncol(dosage) != nsamp) stop("number of dosage columns not equal to number of samples")
      put.var.ncdf(nc, "genotype", dosage, start=c(cnt,1), count=c(nrow(dosage),-1))
      cnt <- cnt + block.size
    }
    close(opfile)

    snps <- as.data.frame(snps, stringsAsFactors=FALSE)
    names(snps) <- c("SNP", "rsID", "position", "alleleA", "alleleB")
  }
  
  if (input.type == "BEAGLE") {
    # .markers file
    if (verbose) message ("Reading SNP file...")
    snps <- read.table(input.files[2], as.is=TRUE, header=FALSE)
    names(snps) <- c("marker", "position", "alleleA", "alleleB")
    if (nrow(snps) != nsnp) stop("SNP number mismatch: ", nsnp, " in genotype file; ", nrow(snps), "in markers file")

    # .dose or .gprobs file
    if (verbose) message ("Reading genotype file...")
    # sample names in header
    opfile <- file(input.files[1], open="r")
    header <- scan(opfile, what=character(), quiet=TRUE, nlines=1)
    header <- header[4:length(header)]
    if (input.dosage) sample.id <- header else sample.id <- header[c(TRUE,FALSE,FALSE)]
    samples <- data.frame("ID"=sample.id, stringsAsFactors=FALSE)
    
    cnt <- 1
    while (length(ss <- scan(opfile, what=character(), quiet=TRUE, nlines=block.size)) > 0) {
      if (verbose) message("Block ", ceiling(cnt/block.size), " of ", ceiling(nsnp/block.size))
      dat <- matrix(ss, ncol=col.cnt, byrow=TRUE)
      if (!allequal(snps[cnt:(cnt+nrow(dat)-1),c(1,3,4)], dat[,1:3])) stop ("markers file does not match genotype file")
      dat <- dat[,4:ncol(dat),drop=FALSE]
      mode(dat) <- "numeric"
      if (input.dosage) {
        # BEAGLE has B allele dosage
        dosage <- 2 - dat
      } else { 
        dosage <- .probToDosage(dat)
      }
      if (ncol(dosage) != nsamp) stop("number of dosage columns not equal to number of samples")
      put.var.ncdf(nc, "genotype", dosage, start=c(cnt,1), count=c(nrow(dosage),-1))
      cnt <- cnt + block.size
    }
    close(opfile)
  }

  if (input.type == "MaCH") {
    # .mlinfo file
    if (verbose) message ("Reading SNP files...")
    snps <- read.table(input.files[2], as.is=TRUE, header=TRUE)
    snps <- snps[,1:3]
    names(snps) <- c("SNP", "alleleA", "alleleB")
    if (nrow(snps) != nsnp) stop("SNP number mismatch: ", nsnp, " in genotype file; ", nrow(snps), "in mlinfo file")

    # file with SNP positions
    snp2 <- read.table(input.files[3], as.is=TRUE, header=TRUE)
    snp2 <- snp2[,c("SNP", "position")]
    if (!setequal(snps$SNP, snp2$SNP)) stop("SNP column in files ", input.files[2], " and ", input.files[3], " do not contain the same values")
    ord <- snps$SNP
    snps <- merge(snps, snp2, sort=FALSE)
    snps <- snps[match(ord, snps$SNP),]

    # empty vector for sample data
    samp.dat <- rep("", nsamp)
    
    # .mldose or .mlprob file
    if (verbose) message ("Reading genotype file...")
    opfile <- file(input.files[1], "r")
    cnt <- 1
    while (length(ss <- scan(opfile, what=character(), quiet=TRUE, nlines=block.size)) > 0) {
      if (verbose) message("Block ", ceiling(cnt/block.size), " of ", ceiling(nsamp/block.size))
      dat <- matrix(ss, ncol=col.cnt, byrow=TRUE)
      samp.dat[cnt:(cnt+nrow(dat)-1)] <- dat[,1]
      dat <- dat[,3:ncol(dat),drop=FALSE]
      mode(dat) <- "numeric"
      if (input.dosage) {
        dosage <- t(dat)
      } else {
        dosage <- t(.probToDosage(dat, BB=FALSE))
      }
      if (nrow(dosage) != nsnp) stop("number of dosage rows not equal to number of SNPs")
      put.var.ncdf(nc, "genotype", dosage, start=c(1,cnt), count=c(-1,ncol(dosage)))
      cnt <- cnt + block.size
    }
    close(opfile)
    
    samp.ids <- strsplit(samp.dat, "->")
    samples <- as.data.frame(matrix(unlist(samp.ids), ncol=2, byrow=TRUE))
    names(samples) <- c("ID_1", "ID_2")
  }
  
  # set up annotation
  samples$scanID <- 1:nrow(samples)
  scanAnnot <- ScanAnnotationDataFrame(samples)
  
  snps$snpID <- 1:nrow(snps)
  # convert chromosome type
  xchr.str <- c(1:22, "X", "Y", "XY", "MT", "M", 23, 24, 25, 26)
  xchr <- as.integer(c(1:22, 23, 25, 24, 26, 26, 23, 24, 25, 26))
  snps$chromosome <- xchr[match(chromosome, xchr.str)]
  snps$position <- as.integer(snps$position)
  snpAnnot <- SnpAnnotationDataFrame(snps)

  # add variable data
  if (verbose) message("Writing annotation...")
  put.var.ncdf(nc, varID, scanAnnot$scanID)
  put.var.ncdf(nc, varpos, snpAnnot$position)
  put.var.ncdf(nc, varchr, snpAnnot$chromosome)
  put.var.ncdf(nc, varA, snpAnnot$alleleA)
  put.var.ncdf(nc, varB, snpAnnot$alleleB)

  # close file
  close.ncdf(nc)

  # save annotation
  save(snpAnnot, file=snp.annot.filename)
  save(scanAnnot, file=scan.annot.filename)
  
  return(invisible(NULL))
}
