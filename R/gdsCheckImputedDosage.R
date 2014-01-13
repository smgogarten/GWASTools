

## genoData needs snpAnnot and scanAnnot attached.
gdsCheckImputedDosage <- function(genoData, snpAnnot, scanAnnot, 
                                  input.files, chromosome,
                                  input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                                  input.dosage=FALSE, block.size=5000,
                                  verbose=TRUE, 
                                  snp.exclude=NULL,
                                  snp.id.start=1,
                                  tolerance=1e-4) {
  
  
  if (missing(snpAnnot)) stop("snp annotation required")
  if (missing(scanAnnot)) stop("scan annotation required")
  
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
  
  # check number of snps and scans
  if (nsnp - length(snp.exclude) != nsnp(genoData)) stop("Number of SNPs does not match.")
  if (nrow(scanAnnot) != nscan(genoData)) stop("Number of scans does not match.")
  

  # read input file(s)
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
    #if (all(is.na(scan.df$sampleID))) {
    #  i_samp <- 1:nrow(samples)
    #  scan.df$sampleID <- samples$ID
    #} else {
    #  i_samp <- match(scan.df$sampleID, samples$ID)
    #}
    
    i_samp <- match(scanAnnot$sampleID, samples$ID)
    
    # empty matrix for SNP data
    snps <- matrix("", nrow=nsnp, ncol=5)
    
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
      
      # add check snp annotation
      stopifnot(allequal(snpAnnot$snp[i_snp : (i_snp+nrow(dat)-1)], dat[, 1]))
      stopifnot(allequal(snpAnnot$rsID[i_snp : (i_snp+nrow(dat)-1)], dat[, 2]))
      stopifnot(allequal(snpAnnot$position[i_snp : (i_snp+nrow(dat)-1)], as.integer(dat[, 3])))
      stopifnot(allequal(snpAnnot$alleleA[i_snp : (i_snp+nrow(dat)-1)], dat[, 4]))
      stopifnot(allequal(snpAnnot$alleleB[i_snp : (i_snp+nrow(dat)-1)], dat[, 5]))
      
      # get dosage info
      dat <- dat[, 6:ncol(dat), drop=FALSE]
      mode(dat) <- "numeric"
      dosage <- .probToDosage(dat)
      if (ncol(dosage) != nsamp) stop("number of dosage columns not equal to number of samples in file")
      
      # subset dosage to match this set of samples, snp subsetting was already taken care of
      dosage <- dosage[, i_samp, drop=FALSE]
      
      # set unphysical dosages to NA
      dosage[dosage < 0 | dosage > 2] <- NA
      
      snp <- c(i_snp, nrow(dosage))
      scan <- c(1, -1)
      
      dosage.geno <- getGenotype(genoData, snp=snp, scan=scan)
      if (!is.matrix(dosage.geno)) dosage.geno <- matrix(dosage.geno, ncol=nscan(genoData), byrow=FALSE)
      
      if (!isTRUE(all.equal(dosage.geno, dosage, tolerance=tolerance))) stop(paste("Dosage not equal between original SNPs", cnt, "and", min(cnt + block.size, nsnp)))
      
      cnt <- cnt + block.size
      i_snp <- i_snp + nrow(dosage)
      
    }
    close(opfile)
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
    
    #if (all(is.na(scan.df$sampleID))) {
    #  i_samp <- 1:nrow(samples)
    #  scan.df$sampleID <- samples$ID
    #} else {
    #  i_samp <- match(scan.df$sampleID, samples$ID)
    #}
    
    i_samp <- match(scanAnnot$sampleID, samples$ID)
    
    # set up snp annotation
    #snp.df[, c("snp", "position", "alleleA", "alleleB")] <- NA
    
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
      
      # add check snp annotation
      stopifnot(allequal(snpAnnot$snp[i_snp : (i_snp+nrow(dat)-1)], dat[, 1]))
      stopifnot(allequal(snpAnnot$alleleA[i_snp : (i_snp+nrow(dat)-1)], dat[, 2]))
      stopifnot(allequal(snpAnnot$alleleB[i_snp : (i_snp+nrow(dat)-1)], dat[, 3]))
      
      
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
      
      # set unphysical dosages to NA
      dosage[dosage < 0 | dosage > 2] <- NA
      
      
      snp <- c(i_snp, nrow(dosage))
      scan <- c(1, -1)
      
      dosage.geno <- getGenotype(genoData, snp=snp, scan=scan)
      if (class(dosage.geno) != "matrix") dosage.geno <- matrix(dosage.geno, ncol=nscan(genoData))
      
      if (!isTRUE(all.equal(dosage.geno, dosage, tolerance=tolerance))) stop(paste("Dosage not equal in original SNPs", cnt, "-", min(nsnp, cnt + block.size)))
      
      cnt <- cnt + block.size

      i_snp <- i_snp + nrow(dosage)
      
    }
    close(opfile)
  }
  
  if (input.type == "MaCH") {
    # .mlinfo file
    if (verbose) message ("Reading SNP files...")
    snps <- read.table(input.files[2], as.is=TRUE, header=TRUE)
    snps <- snps[,1:3]
    names(snps) <- c("snp", "alleleA", "alleleB")
    if (nrow(snps) != nsnp) stop("SNP number mismatch: ", nsnp, " in genotype file; ", nrow(snps), "in mlinfo file")
    
    # file with SNP positions
    snp2 <- read.table(input.files[3], as.is=TRUE, header=TRUE)
    snp2 <- snp2[,c("SNP", "position")]
    names(snp2) <- c("snp", "position")
    if (!setequal(snps$snp, snp2$snp)) stop("SNP column in files ", input.files[2], " and ", input.files[3], " do not contain the same values")

    # get snps to exclude
    i_snp <- !((1:nsnp) %in% snp.exclude)
    
    # check snp annotation
    stopifnot(allequal(snpAnnot$snp, snps$snp[i_snp]))
    stopifnot(allequal(snpAnnot$alleleA, snps$alleleA[i_snp]))
    stopifnot(allequal(snpAnnot$alleleB, snps$alleleB[i_snp]))
    stopifnot(allequal(snpAnnot$position, snp2$position[i_snp]))
    
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
      # i_snp matches the ordering of snps in the dosage array to the ordering of snps in the gds file
      for (i_dos in 1:length(samp.block)) {
        
        # this sample
        samp <- samp.block[i_dos]
        
        i_samp <- which(scanAnnot$sampleID %in% samp)
        if (length(i_samp) < 1) next
        
        # get dosage for that sample and reorder to match gds snp ordering.
        dos.samp <- dosage[i_snp, i_dos]
        
        # set unphysical dosages to NA
        dos.samp[dos.samp < 0 | dos.samp > 2] <- NA
        
        snp <- c(1, -1)
        scan <- c(i_samp, 1)
        
        dosage.geno <- getGenotype(genoData, snp=snp, scan=scan)

        if (!isTRUE(all.equal(dos.samp, dosage.geno, tolerance=tolerance))) stop(paste("Dosage not equal between original sample", cnt, "and", min(cnt+block.size, nsamp)))
        
        cnt <- cnt+1
      }
    }
    close(opfile)
  }
  
  scanID.gds <- getScanID(genoData)
  
  # check that all added=FALSE samples are NA
  i_zero <- scanAnnot$scanID[!scanAnnot$added]
  
  for (i_samp in i_zero) {
    # match gds file
    i_gds <- match(i_samp, scanID.gds)
    geno.gds <- getGenotype(genoData, snp=c(1,-1), scan=c(i_gds, 1))
    if (!(all(is.na(geno.gds)))) stop(paste("genotypes are not NA for added=FALSE sample", i_samp))
  }
  
  
  return(TRUE)
}