.probToDosage <- function(probs, BB=TRUE, prob.miss.val=NULL) {
  if (BB & ncol(probs) %% 3 != 0) stop("invalid probability file - there are not 3 columns per row")
  if (!BB & ncol(probs) %% 2 != 0) stop("invalid probability file - there are not 2 columns per row")

  if (BB) {
    AAprob <- probs[,c(TRUE,FALSE,FALSE),drop=FALSE]
    ABprob <- probs[,c(FALSE,TRUE,FALSE),drop=FALSE]
  } else {
    AAprob <- probs[,c(TRUE,FALSE),drop=FALSE]
    ABprob <- probs[,c(FALSE,TRUE),drop=FALSE]
  }

  # calculate A allele dosage
  dosage <- 2*AAprob + ABprob
  
  return(dosage)
}


gdsImputedDosage <- function(input.files, gds.filename, chromosome,
                             input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                             input.dosage=FALSE, block.size=5000,
                             snp.annot.filename="dosage.snp.RData",
                             scan.annot.filename="dosage.scan.RData",
                             verbose=TRUE, zipflag="ZIP.max", genotypeDim="snp,scan",
                             scan.df=NULL,
                             snp.exclude=NULL,
                             snp.id.start=1) {
  
  # arguments: input type (impute2, beagle, mach), probs or dosages
  input.type <- match.arg(input.type)
  
  # check zipflag
  if (!(zipflag %in% c("", "ZIP", "ZIP.fast", "ZIP.default", "ZIP.max"))) stop("zipflag must be one of ZIP, ZIP.fast, ZIP.default, ZIP.max")
  
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

  
  # assign snp IDs
  if (!is.null(snp.exclude)){
    # check type
    if (!is.numeric(snp.exclude)) stop("snp.exclude must be list of integers to include")
    snpID <- 1:(nsnp - length(snp.exclude))
  } else {
    if (verbose) message("Including all SNPs.")
    snpID <- 1:nsnp
    snp.exclude <- numeric()
    #snp.names <- NA
  }
  snp.df <- data.frame(snpID=as.integer(snpID + (snp.id.start-1)))
  
  
  # scan.df - input checking or creation
  if (!is.null(scan.df) & !("sampleID" %in% names(scan.df))) {
    stop("sampleID required in scan.df.")
    if (!("scanID" %in% names(scan.df))) scan.df$scanID <- 1:nrow(scan.df)
  } else if (!is.null(scan.df) & !("scanID" %in% names(scan.df))) {
    if (verbose) message("scanID mapping not given. Assigning scanIDs automatically.")
    scan.df$scanID <- 1:nsamp
  } else if (is.null(scan.df)) {
    if (verbose) message("scan.df not given. Assigning scanIDs automatically.")
    scan.df <- data.frame(scanID=1:nsamp, stringsAsFactors=FALSE)
    scan.df$sampleID <- NA
  }
  # order by scanID -- do we want this?
  #scan.df$scanID <- as.integer(scan.df$scanID) # convert to integer
  #scan.df <- scan.df[order(scan.df$scanID), , drop=FALSE]
  scan.df$added <- FALSE # variable if added to the gds file successfully
  #if (nrow(scan.df) != nsamp) stop("Number of samples in scan.df is different than number of samples in file.")
  
  # add some checking so we can use our own scan ids (and our own snp ids?)
  #snpID <- 1:nsnp
  #scanID <- 1:nsamp
  
  # create GDS
  gfile <- createfn.gds(gds.filename)
  
  add.gdsn(gfile, "snp.id", snp.df$snpID, compress=zipflag, closezip=TRUE)
  add.gdsn(gfile, "sample.id", scan.df$scanID, compress=zipflag, closezip=TRUE)
  
  n <- add.gdsn(gfile, name="description")
  put.attr.gdsn(n, "FileFormat", "IMPUTED_DOSAGE")
  
  geno.valdim <- switch(genotypeDim,
                        "snp,scan"=c(length(snp.df$snpID), length(scan.df$scanID)),
                        "scan,snp"=c(length(scan.df$scanID), length(snp.df$snpID)))
  gGeno <- add.gdsn(gfile, "genotype", valdim=geno.valdim, storage="float32")
  
  geno.order <- switch(genotypeDim,
                       "snp,scan"="snp.order",
                       "scan,snp"="sample.order")
  put.attr.gdsn(gGeno, geno.order)
  miss.val <- -1
  put.attr.gdsn(gGeno, "missing.value", miss.val)
  
  
  # read input file(s)
  # valid for GDS and NCDF mostly.
  if (input.type == "IMPUTE2") {
    
    if (input.dosage) stop("input.dosage=TRUE not valid for input.type=IMPUTE2")
    # .samples file
    if (verbose) message ("Reading sample file...")
    samp.header <- scan(input.files[2], what=character(), nlines=1, quiet=TRUE)
    samples <- read.table(input.files[2], as.is=TRUE, header=FALSE, skip=2)
    names(samples) <- samp.header
    if (nrow(samples) != nsamp) stop("Sample number mismatch: ", nsamp, " in genotype file; ", nrow(samples), "in samples file")

    # check for unique ids
    samples$ID <- paste(samples$ID_1, samples$ID_2)
    if (any(duplicated(samples$ID))) stop("Sample ID is duplicated in input sample file")
    
    # select samples
    if (all(is.na(scan.df$sampleID))) {
      i_samp <- 1:nrow(samples)
      scan.df$sampleID <- samples$ID
    } else {
      i_samp <- match(scan.df$sampleID, samples$ID)
    }

    # empty matrix for SNP data
    snps <- matrix("", nrow=nsnp, ncol=5)
    
    # set up snp df columns
    snp.df[, c("snp", "rsID", "position", "alleleA", "alleleB")] <- NA
    
    # .gens file
    if (verbose) message ("Reading genotype file...")
    opfile <- file(input.files[1], "r")
    cnt <- 1 # tracks where we are in the dosage file, snp dimension
    i_snp <- 1 # tracks where we are in the gds file, snp dimension
    while (length(ss <- scan(opfile, what=character(), quiet=TRUE, nlines=block.size)) > 0) {
      if (verbose) message("Block ", ceiling(cnt/block.size), " of ", ceiling(nsnp/block.size))
      dat <- matrix(ss, ncol=col.cnt, byrow=TRUE) # snps are the first column
      
      # selected snps
      dat <- dat[!((1:nrow(dat)+cnt-1) %in% snp.exclude), , drop=FALSE]
      # check that all are selected
      if (nrow(dat) < 1){
        cnt <- cnt + block.size
        next
      }
      
      # add snp info to annotation
      snp.df$snp[i_snp : (i_snp+nrow(dat)-1)] <- dat[, 1]
      snp.df$rsID[i_snp : (i_snp+nrow(dat)-1)] <- dat[, 2]
      snp.df$position[i_snp : (i_snp+nrow(dat)-1)] <- as.integer(dat[, 3])
      snp.df$alleleA[i_snp : (i_snp+nrow(dat)-1)] <- dat[, 4]
      snp.df$alleleB[i_snp : (i_snp+nrow(dat)-1)] <- dat[, 5]
      
      # get dosage info
      dat <- dat[, 6:ncol(dat), drop=FALSE]
      mode(dat) <- "numeric"
      dosage <- .probToDosage(dat)
      if (ncol(dosage) != nsamp) stop("number of dosage columns not equal to number of samples in file")
      
      # subset dosage to match this set of samples, snp subsetting was already taken care of
      dosage <- dosage[, i_samp, drop=FALSE]
      if(genotypeDim == "snp,scan") {
        start <- c(i_snp, 1)
        count <- c(nrow(dosage), ncol(dosage))
      } else if (genotypeDim == "scan,snp") {
        start <- c(1, i_snp)
        count <- c(ncol(dosage), nrow(dosage))
        dosage <- t(dosage)
      }
      
      # check for missing values
      dosage[dosage < 0 | dosage > 2] <- miss.val

      write.gdsn(gGeno, dosage, start=start, count=count)
      cnt <- cnt + block.size
      if (genotypeDim == "snp,scan"){
        i_snp <- i_snp + nrow(dosage)
      } else if (genotypeDim == "scan,snp"){
        i_snp <- i_snp + ncol(dosage)
      }
    
    }
    close(opfile)

    scan.df$added[scan.df$sampleID %in% samples$ID] <- TRUE

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
    
    if (all(is.na(scan.df$sampleID))) {
      i_samp <- 1:nrow(samples)
      scan.df$sampleID <- samples$ID
    } else {
      i_samp <- match(scan.df$sampleID, samples$ID)
    }
    
    # set up snp annotation
    snp.df[, c("snp", "position", "alleleA", "alleleB")] <- NA
    
    cnt <- 1 # tracks location in the impute file
    i_snp <- 1 # tracks location in the gds file
    while (length(ss <- scan(opfile, what=character(), quiet=TRUE, nlines=block.size)) > 0) {
      
      if (verbose) message("Block ", ceiling(cnt/block.size), " of ", ceiling(nsnp/block.size))
      
      dat <- matrix(ss, ncol=col.cnt, byrow=TRUE)
      
      if (!allequal(snps[cnt:(cnt+nrow(dat)-1),c(1,3,4)], dat[,1:3])) stop ("markers file does not match genotype file")
      
      # selected snps
      dat <- dat[!((1:nrow(dat)+cnt-1) %in% snp.exclude), , drop=FALSE]
      # check that all are selected
      if (nrow(dat) < 1){
        cnt <- cnt + block.size
        next
      }
      
      # add snp info to annotation
      snp.df$snp[i_snp : (i_snp + nrow(dat) - 1)] <- dat[, 1]
      snp.df$alleleA[i_snp : (i_snp + nrow(dat) - 1)] <- dat[, 2]
      snp.df$alleleB[i_snp : (i_snp + nrow(dat) - 1)] <- dat[, 3]
      
      dat <- dat[, 4:ncol(dat), drop=FALSE]
      mode(dat) <- "numeric"
      if (input.dosage) {
        # BEAGLE has B allele dosage
        dosage <- 2 - dat
      } else { 
        dosage <- .probToDosage(dat)
      }
      if (ncol(dosage) != nsamp) stop("number of dosage columns not equal to number of samples")
      dosage <- dosage[, i_samp, drop=FALSE]
      if(genotypeDim == "snp,scan") {
        start <- c(i_snp, 1)
        count <- c(nrow(dosage), -1)
      } else if (genotypeDim == "scan,snp") {
        start <- c(1, i_snp)
        count <- c(-1, nrow(dosage))
        dosage <- t(dosage)
      }

      # check for missing values
      dosage[dosage < 0 | dosage > 2] <- miss.val
      
      write.gdsn(gGeno, dosage, start=start, count=count)
      cnt <- cnt + block.size
      if (genotypeDim == "snp,scan"){
        i_snp <- i_snp + nrow(dosage)
      } else if (genotypeDim == "scan,snp"){
        i_snp <- i_snp + ncol(dosage)
      }
      
    }
    close(opfile)
    
    # keep track of added samples
    scan.df$added[scan.df$sampleID %in% samples$ID] <- TRUE

    # add position
    snp.df$position <- snps$position[!((1:nrow(snps)) %in% snp.exclude)]
  }

  if (input.type == "MaCH") {
    # .mlinfo file
    if (verbose) message ("Reading SNP files...")
    snps <- read.table(input.files[2], as.is=TRUE, header=TRUE)
    snps <- snps[,1:3]
    names(snps) <- c("snp", "alleleA", "alleleB")
    if (nrow(snps) != nsnp) stop("SNP number mismatch: ", nsnp, " in genotype file; ", nrow(snps), "in mlinfo file")
        
    # get snps to exclude
    i_snp <- !((1:nsnp) %in% snp.exclude)
    
    snp.df[, c("snp", "alleleA", "alleleB")] <- NA
    snp.df$snp <- snps$snp[i_snp]
    snp.df$alleleA <- snps$alleleA[i_snp]
    snp.df$alleleB <- snps$alleleB[i_snp]
    
    # file with SNP positions
    snp2 <- read.table(input.files[3], as.is=TRUE, header=TRUE)
    snp2 <- snp2[,c("SNP", "position")]
    names(snp2) <- c("snp", "position")
    if (!setequal(snps$snp, snp2$snp)) stop("SNP column in files ", input.files[2], " and ", input.files[3], " do not contain the same values")
    #ord <- snps$SNP
    #snp.df <- merge(snp.df, snp2, sort=FALSE, all.x=TRUE)
    #snp.df <- snp.df[order(snp.df$snpID),]
    snp.df$position <- snp2$position[i_snp]

    
    # .mldose or .mlprob file
    if (verbose) message ("Reading genotype file...")
    opfile <- file(input.files[1], "r")
    cnt <- 1
    while (length(ss <- scan(opfile, what=character(), quiet=TRUE, nlines=block.size)) > 0) {
      if (verbose) message("Block ", ceiling(cnt/block.size), " of ", ceiling(nsamp/block.size))
      dat <- matrix(ss, ncol=col.cnt, byrow=TRUE)
      #samp.dat[cnt:(cnt+nrow(dat)-1)] <- dat[,1]
      samp.block <- dat[,1]
      dat <- dat[,3:ncol(dat),drop=FALSE]
      mode(dat) <- "numeric"
      if (input.dosage) {
        dosage <- t(dat)
      } else {
        dosage <- t(.probToDosage(dat, BB=FALSE))
      }
      if (nrow(dosage) != nsnp) stop("number of dosage rows not equal to number of SNPs")
            
      # loop over samples to add them. Lots of indices here:
      # i_dos tracks the location in the dosage array
      # i_samp finds the location in the gds file/scan.df
      # i_snp matches the ordering of snps in the dosage array to the ordering of snps in the gdsf ile
      for (i_dos in 1:length(samp.block)) {
        
        # this sample
        samp <- samp.block[i_dos]
        
        # check if sampleID has been given
        if (all(!is.na(scan.df$sampleID)) & (samp %in% scan.df$sampleID)) {
          # add the genotype
          i_samp <- which(scan.df$sampleID %in% samp)
        } else if (any(is.na(scan.df$sampleID)) & !(samp %in% scan.df$sampleID)) {
          # take the first missing sampleID
          i_samp <- min(which(is.na(scan.df$sampleID)))
          scan.df$sampleID[i_samp] <- samp
        } else {
          # something failed?
          next
        }
        # get dosage for that sample and reorder to match gds snp ordering.
        dos.samp <- dosage[i_snp, i_dos]
        
        # check for missing values
        dos.samp[dos.samp < 0 | dos.samp > 2] <- miss.val
        
        if (genotypeDim == "snp,scan"){
          start <- c(1, i_samp)
          count <- c(-1, 1)
        } else if (genotypeDim == "scan,snp") {
          start <- c(i_samp, 1)
          count <- c(1, -1)
        }
        # write all snps for that sample
        write.gdsn(gGeno, dos.samp, start=start, count=count)
        scan.df$added[i_samp] <- TRUE
        cnt <- cnt+1
      }
    }
#       if(genotypeDim == "snp,scan") {
#         start <- c(1, cnt)
#         count <- c(-1, ncol(dosage))
#       } else if (genotypeDim == "scan,snp") {
#         start <- c(cnt, 1)
#         count <- c(ncol(dosage), -1)
#         dosage <- t(dosage)
#       }
# 
#       write.gdsn(gGeno, dosage, start=start, count=count)
#       cnt <- cnt + block.size

    close(opfile)
    

    #tmp <- as.data.frame(matrix(unlist(strsplit(scan.df$sampleID, "->")), ncol=2, byrow=TRUE))
    #names(tmp) <- c("ID_1", "ID_2")
    #scan.df$ID_1 <- tmp$ID_1
    #scan.df$ID_2 <- tmp$ID_2
  }

  # zero out samples that haven't been added
  i_zero <- which(!scan.df$added)
  for (i_samp in i_zero){
    if (genotypeDim == "snp,scan"){
      start <- c(1, i_samp)
      count <- c(-1, 1)
    } else if (genotypeDim == "scan,snp") {
      start <- c(i_samp, 1)
      count <- c(1, -1)
    }
    dos.samp <- rep(-1, nrow(snp.df))
    # write all snps for that sample
    write.gdsn(gGeno, dos.samp, start=start, count=count)
  }
  
  # set up annotation
  #samples$scanID <- 1:nrow(samples) ## this will need to change.
  if ("sex" %in% names(samples)) names(samples)[names(samples) %in% "sex"] <- "sex.sample"
  scanAnnot <- ScanAnnotationDataFrame(scan.df)
  
  #snps$snpID <- 1:nrow(snps)
  # convert chromosome type
  xchr.str <- c(1:22, "X", "Y", "XY", "MT", "M", 23, 24, 25, 26)
  xchr <- as.integer(c(1:22, 23, 25, 24, 26, 26, 23, 24, 25, 26))
  snp.df$chromosome <- xchr[match(chromosome, xchr.str)]
  snp.df$position <- as.integer(snp.df$position)
  snpAnnot <- SnpAnnotationDataFrame(snp.df)

  # add variable data
  if (verbose) message("Writing annotation...")
  # add "sample.id"
  #add.gdsn(gfile, "sample.id", scanAnnot$scanID, compress=zipflag, closezip=TRUE)  
  # add "snp.id"
  #add.gdsn(gfile, "snp.id", snpAnnot$snpID, compress=zipflag, closezip=TRUE)
  # add "snp.chromosome"
  add.gdsn(gfile, "snp.chromosome", snpAnnot$chromosome, compress=zipflag, closezip=TRUE)
  # add "snp.position"
  add.gdsn(gfile, "snp.position", snpAnnot$position, compress=zipflag, closezip=TRUE)
  # add alleles
  add.gdsn(gfile, "snp.allele", paste(snpAnnot$alleleA, snpAnnot$alleleB, sep="/"), compress=zipflag, closezip=TRUE)

  
  # close file
  #close.ncdf(nc)
  closefn.gds(gfile)
  
  # clean up gds file
  cleanup.gds(gds.filename, verbose=verbose)
  
  # save annotation
  save(snpAnnot, file=snp.annot.filename)
  save(scanAnnot, file=scan.annot.filename)
  
  return(invisible(NULL))
}

