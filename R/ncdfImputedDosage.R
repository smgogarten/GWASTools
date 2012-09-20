.probToDosage <- function(probs, BB=TRUE) {
  if (BB & ncol(probs) %% 3 != 0) stop("invalid probability file - there are not 3 columns per row")
  if (!BB & ncol(probs) %% 2 != 0) stop("invalid probability file - there are not 2 columns per row")

  if (BB) {
    AAprob <- probs[,c(TRUE,FALSE,FALSE)]
    ABprob <- probs[,c(FALSE,TRUE,FALSE)]
  } else {
    AAprob <- probs[,c(TRUE,FALSE)]
    ABprob <- probs[,c(FALSE,TRUE)]
  }

  # calculate A allele dosage
  dosage <- 2*AAprob + ABprob

  return(dosage)
}

ncdfImputedDosage <- function(input.files, ncdf.filename, chromosome,
                              input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                              input.dosage=FALSE, block.size=5000,
                              snp.annot.filename="dosage.snp.RData",
                              scan.annot.filename="dosage.scan.RData") {
  # arguments: input type (impute2, beagle, mach), probs or dosages

  # loop, reading rows in blocks for large files
  # read input file(s)
  input.type <- match.arg(input.type)
  if (input.type == "IMPUTE2") {
    if (input.dosage) stop("input.dosage=TRUE not valid for input.type=IMPUTE2")
    # .samples file
    samp.header <- scan(input.files[2], what=character(), nlines=1, quiet=TRUE)
    samples <- read.table(input.files[2], as.is=TRUE, header=FALSE, skip=2)
    names(samples) <- samp.header
    
    # .gens file
    dat <- read.table(input.files[1], as.is=TRUE, header=FALSE)
    snps <- dat[,1:5]
    names(snps) <- c("SNP", "rsID", "position", "alleleA", "alleleB")
    dat <- dat[,6:ncol(dat)]
    dosage <- .probToDosage(as.matrix(dat))
  }
  
  if (input.type == "BEAGLE") {
    #opfile <- file(input.files[1], open="r")
    #header <- scan(text=readLines(opfile, n=1), what=character(), quiet=TRUE)
    #close(opfile)
    
    # .markers file
    snps <- read.table(input.files[2], as.is=TRUE, header=FALSE)
    names(snps) <- c("marker", "position", "alleleA", "alleleB")

    # .dose or .gprobs file
    dat <- read.table(input.files[1], as.is=TRUE, header=TRUE)
    stopifnot(allequal(snps[,c(1,3,4)], dat[,1:3,]))
    dat <- dat[,4:ncol(dat)]

    if (input.dosage) {
      # .dose file
      sample.id <- names(dat)
      # BEAGLE has B allele dosage
      dosage <- 2 - as.matrix(dat)
    } else {      
      # .grobs file
      sample.id <- names(dat)[c(TRUE,FALSE,FALSE)]
      dosage <- .probToDosage(as.matrix(dat))
    }
    samples <- data.frame("ID"=sample.id, stringsAsFactors=FALSE)
  }

  if (input.type == "MaCH") {
    # .mlinfo file
    snps <- read.table(input.files[2], as.is=TRUE, header=TRUE)
    snps <- snps[,1:3]
    names(snps) <- c("SNP", "alleleA", "alleleB")

    # file with SNP positions
    snp2 <- read.table(input.files[3], as.is=TRUE, header=TRUE)
    snp2 <- snp2[,c("SNP", "position")]
    if (!setequal(snps$SNP, snp2$SNP)) stop("SNP column in files ", input.files[2], " and ", input.files[3], " do not contain the same values")
    ord <- snps$SNP
    snps <- merge(snps, snp2, sort=FALSE)
    snps <- snps[match(ord, snps$SNP),]

    # .mldose or .mlprob file
    dat <- read.table(input.files[1], as.is=TRUE, header=FALSE)
    samp.ids <- strsplit(dat[,1], "->")
    samples <- as.data.frame(matrix(unlist(samp.ids), ncol=2, byrow=TRUE))
    names(samples) <- c("ID_1", "ID_2")
    dat <- dat[,3:ncol(dat)]

    if (input.dosage) {
      # .mldose file
      dosage <- t(as.matrix(dat))
    } else {
      dosage <- t(.probToDosage(as.matrix(dat), BB=FALSE))
    }
  }
  
#  while (length(s <- readLines(opfile, n=block.size)) > 0) {
#    for (i in 1:length(s)){ 
#      ss <- scan(text=s[i], what=character(), quiet=TRUE)      
#    }
#  }
  
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

  # create NetCDF
  # define dimensions
  snpdim <- dim.def.ncdf("snp", "count", snpAnnot$snpID)
  sampledim <- dim.def.ncdf("sample", "count", 1:nrow(scanAnnot), unlim=TRUE)
  chardim <- dim.def.ncdf("nchar", "", 1)

  # define variables
  varID <- var.def.ncdf("sampleID", "id", dim=sampledim, missval=0, prec="integer")
  varpos <- var.def.ncdf("position", "bases", dim=snpdim, missval=-1, prec="integer")
  varchr <- var.def.ncdf("chromosome", "id", dim=snpdim, missval=-1, prec="integer")
  varA <- var.def.ncdf("alleleA", "allele", dim=list(chardim,snpdim), missval="0", prec="char")
  varB <- var.def.ncdf("alleleB", "allele", dim=list(chardim,snpdim), missval="0", prec="char")
  vargeno <- var.def.ncdf("genotype", "A_allele_dosage", dim=list(snpdim,sampledim), missval=-1, prec="single")

  # create the NetCDF file
  ncfile <- create.ncdf(ncdf.filename, list(varID, varpos, varchr, varA, varB, vargeno))
        
  # add variable data
#  if (verbose) message(date(), "\t\twriting position and chromosome ...\n")
  put.var.ncdf(ncfile, varpos, snpAnnot$position)
  put.var.ncdf(ncfile, varchr, snpAnnot$chromosome)
  put.var.ncdf(ncfile, varA, snpAnnot$alleleA)
  put.var.ncdf(ncfile, varB, snpAnnot$alleleB)

  # add genotype data
#  if (verbose) message(date(), "\t\twriting genotypes ...\n")
  put.var.ncdf(ncfile, "genotype", dosage)

  # add sample id
  put.var.ncdf(ncfile, varID, scanAnnot$scanID)

  # close file
  close.ncdf(ncfile)

  # save annotation
  save(snpAnnot, file=snp.annot.filename)
  save(scanAnnot, file=scan.annot.filename)
  
  return(invisible(NULL))
}
