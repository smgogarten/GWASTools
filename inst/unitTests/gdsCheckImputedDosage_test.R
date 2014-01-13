
test_check_beagle <- function() {
  probfile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.gprobs",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.dose",
                          package="GWASdata")
  markfile <- system.file("extdata", "imputation", "BEAGLE", "hapmap.markers",
                          package="GWASdata")
  
  # subsets of data to read in
  header <- read.table(dosefile, header=FALSE, stringsAsFactors=FALSE, nrow=1)
  
  markers <- read.table(markfile, header=FALSE, stringsAsFactors=FALSE)
  
  
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  nsnp <- nrow(markers)
  nscan <- ncol(header)-3
  
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(1:nscan, 1)
  
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 40, 99)
  genoDim <- "snp,scan"
    for (b in blocks) {
      for (i in 1:2) {
        gdsImputedDosage(input.files=c(files[i], markfile), gds.filename=gdsfile, chromosome=22,
                         input.type="BEAGLE", input.dosage=inputs[i], block.size=b,
                         snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim,
                         verbose=FALSE)
        
        gds <- GdsGenotypeReader(gdsfile)
        scanAnnot <- getobj(scanfile)
        snpAnnot <- getobj(snpfile)
        genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

        gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                              input.files=c(files[i], markfile), chromosome,
                              input.type="BEAGLE", 
                              input.dosage=inputs[i], block.size=b)
        
        geno.orig <- getGenotype(genoData) # for checking later
        close(genoData)
        
        # change a snp
        gds <- openfn.gds(gdsfile, readonly=FALSE)
        orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
        val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
        write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
        closefn.gds(gds)
        
        # check exception
        gds <- GdsGenotypeReader(gdsfile)
        genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
        geno <- getGenotype(genoData)
        
        checkException(checkEquals(geno.orig, geno)) # make sure they actualy are different
        checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                             input.files=c(files[i], markfile), chromosome,
                                             input.type="BEAGLE", 
                                             input.dosage=inputs[i], block.size=b))
        
        close(genoData)
      }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}


test_check_beagle_missing <- function() {
  probfile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.gprobs",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.dose",
                          package="GWASdata")
  markfile <- system.file("extdata", "imputation", "BEAGLE", "hapmap.markers",
                          package="GWASdata")
  
  # subsets of data to read in
  header <- read.table(dosefile, header=FALSE, stringsAsFactors=FALSE, nrow=1)
  
  markers <- read.table(markfile, header=FALSE, stringsAsFactors=FALSE)
  
  newprobfile <- tempfile()  
  prob <- read.table(probfile, header=TRUE, stringsAsFactors=FALSE)
  prob[1, 4:6] <- -1
  write.table(prob, file=newprobfile, row.names=FALSE, col.names=TRUE)
  x <- read.table(newprobfile, header=TRUE, stringsAsFactors=FALSE)
  checkEquals(prob, x)
  
  newdosefile <- tempfile()
  dose <- read.table(dosefile, header=TRUE, stringsAsFactors=FALSE)
  dose[1, 4] <- -1
  write.table(dose, file=newdosefile, row.names=FALSE, col.names=TRUE)
  x <- read.table(newdosefile, header=TRUE, stringsAsFactors=FALSE)
  checkEquals(dose, x)
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  nsnp <- nrow(markers)
  nscan <- ncol(header)-3
  
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(1:nscan, 1)
  
  
  files <- c(newprobfile, newdosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 40, 99)
  genoDim <- "snp,scan"
  for (b in blocks) {
    for (i in 1:2) {
      gdsImputedDosage(input.files=c(files[i], markfile), gds.filename=gdsfile, chromosome=22,
                       input.type="BEAGLE", input.dosage=inputs[i], block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim,
                       verbose=FALSE)
      
      gds <- GdsGenotypeReader(gdsfile)
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      
      gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                            input.files=c(files[i], markfile), chromosome,
                            input.type="BEAGLE", 
                            input.dosage=inputs[i], block.size=b)
      
      geno.orig <- getGenotype(genoData) # for checking later
      close(genoData)
      
      # change a snp
      gds <- openfn.gds(gdsfile, readonly=FALSE)
      orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
      val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
      write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
      closefn.gds(gds)
      
      # check exception
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      geno <- getGenotype(genoData)
      
      checkException(checkEquals(geno.orig, geno)) # make sure they actualy are different
      checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                           input.files=c(files[i], markfile), chromosome,
                                           input.type="BEAGLE", 
                                           input.dosage=inputs[i], block.size=b))
      
      close(genoData)
    }
  }
  
  unlink(c(gdsfile, snpfile, scanfile))
}


# tests adding only a subset of samples/snps
test_check_beagle_subset <- function() {
  probfile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.gprobs",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.dose",
                          package="GWASdata")
  markfile <- system.file("extdata", "imputation", "BEAGLE", "hapmap.markers",
                          package="GWASdata")
  
  # subsets of data to read in
  header <- read.table(dosefile, header=FALSE, stringsAsFactors=FALSE, nrow=1)
  # reverse the samples and insert a missing sample between
  scan.df <- data.frame(sampleID=c(header[1,5], "missing", header[1,4]), scanID=c(200, 205, 208))
  markers <- read.table(markfile, header=FALSE, stringsAsFactors=FALSE)
  i_snp_rm <- sample(1:nrow(markers), 5)
  #snp.names <- markers$V1[i_snp_rm] # remove 5 random SNPs
  nsnp <- nrow(markers) - length(i_snp_rm)
  
  snp.id.start <- 100
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  # random indices to change
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(which(scan.df$sampleID != "missing"), 1)
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 40, 99, 1)
  genoDim <- "snp,scan"
  for (b in blocks) {
    for (i in 1:2) {
      
      # test reading snp.names and scan.df
      gdsImputedDosage(input.files=c(files[i], markfile), gds.filename=gdsfile, chromosome=22,
                       input.type="BEAGLE", input.dosage=inputs[i], block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                       scan.df=scan.df, snp.exclude=i_snp_rm, genotypeDim=genoDim,
                       snp.id.start=snp.id.start, verbose=FALSE)
      
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      
      gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                            input.files=c(files[i], markfile), chromosome,
                            input.type="BEAGLE", 
                            input.dosage=inputs[i], block.size=b,
                            snp.exclude=i_snp_rm)
      
      # change it in a random place
      close(genoData)
      
      gds <- openfn.gds(gdsfile, readonly=FALSE)
      orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
      val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
      write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
      closefn.gds(gds)
      
      # check exception
      
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                            input.files=c(files[i], markfile), chromosome,
                            input.type="BEAGLE", 
                            input.dosage=inputs[i], block.size=b,
                            snp.exclude=i_snp_rm))
      
      close(genoData)
      
      
      # change an added=FALSE sample
      tmp <- which(!scanAnnot$added)
      i_scan2 <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
      i_snp2 <- sample(1:nsnp, 1)
      gds <- openfn.gds(gdsfile, readonly=FALSE)
      # put the original value back
      write.gdsn(index.gdsn(gds, "genotype"), orig, start=c(i_snp,i_scan), count=c(1,1))
      # write 
      write.gdsn(index.gdsn(gds, "genotype"), 1, start=c(i_snp2, i_scan2), count=c(1,1))
      closefn.gds(gds)
      
      # check exception
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                           input.files=c(files[i], markfile), chromosome,
                                           input.type="BEAGLE", 
                                           input.dosage=inputs[i], block.size=b,
                                           snp.exclude=i_snp_rm))
      

      close(genoData)
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}



test_check_mach <- function() {
  probfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlprob",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mldose",
                          package="GWASdata")
  markfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlinfo",
                          package="GWASdata")
  posfile <- system.file("extdata", "imputation", "MaCH", "mach1.snp.position",
                         package="GWASdata")
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  # set up scan.df for subsetting later.
  dosages <- read.table(dosefile, header=FALSE, stringsAsFactors=FALSE)
  
  markers <- read.table(markfile, header=TRUE, stringsAsFactors=FALSE)
  
  nsnp <- nrow(markers)
  nscan <- nrow(dosages)
  
  
  # random snp and scan to remove
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(1:nscan, 1)
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 40, 99)
  genoDim <- "snp,scan"
    for (b in blocks) {
      for (i in 1:2) {
        gdsImputedDosage(input.files=c(files[i], markfile, posfile), gds.filename=gdsfile, chromosome=22,
                         input.type="MaCH", input.dosage=inputs[i], block.size=b,
                         snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                         genotypeDim=genoDim, verbose=FALSE)
        
        gds <- GdsGenotypeReader(gdsfile)
        scanAnnot <- getobj(scanfile)
        snpAnnot <- getobj(snpfile)
        genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
        
        gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                              input.files=c(files[i], markfile, posfile), chromosome,
                              input.type="MaCH", 
                              input.dosage=inputs[i], block.size=b)
        
        geno.orig <- getGenotype(genoData)
        
        close(genoData)
        
        # change a snp
        gds <- openfn.gds(gdsfile, readonly=FALSE)
        orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
        val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
        write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
        closefn.gds(gds)
        
        # check exception
        gds <- GdsGenotypeReader(gdsfile)
        genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
        
        geno <- getGenotype(genoData)
        
        checkException(checkEquals(geno, geno.orig))
        checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                             input.files=c(files[i], markfile, posfile), chromosome,
                                             input.type="MaCH", 
                                             input.dosage=inputs[i], block.size=b))
        
        
        close(genoData)
      }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}


test_check_mach_missing <- function() {

  probfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlprob",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mldose",
                          package="GWASdata")
  markfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlinfo",
                          package="GWASdata")
  posfile <- system.file("extdata", "imputation", "MaCH", "mach1.snp.position",
                         package="GWASdata")
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  newprobfile <- tempfile()
  newdosefile <- tempfile()
  
  dosages <- read.table(dosefile, header=FALSE, stringsAsFactors=FALSE)
  dosages[1, 3] <- -1
  write.table(dosages, file=newdosefile, row.names=FALSE, col.names=FALSE)
  x <- read.table(newdosefile, stringsAsFactors=FALSE, header=FALSE)
  checkEquals(x, dosages)
  
  # mach probs are AA, AB not AA, AB, BB
  prob <- read.table(probfile, header=FALSE, stringsAsFactors=FALSE)
  prob[1, 3:4]  <- -1
  write.table(prob, file=newprobfile, row.names=FALSE, col.names=FALSE)
  x <- read.table(newprobfile, stringsAsFactors=FALSE, header=FALSE)
  checkEquals(x, prob)
  
  markers <- read.table(markfile, header=TRUE, stringsAsFactors=FALSE)
  
  nsnp <- nrow(markers)
  nscan <- nrow(dosages)
  
  
  # random snp and scan to remove
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(1:nscan, 1)
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 40, 99)
  genoDim <- "snp,scan"
  for (b in blocks) {
    for (i in 1:2) {
      gdsImputedDosage(input.files=c(files[i], markfile, posfile), gds.filename=gdsfile, chromosome=22,
                       input.type="MaCH", input.dosage=inputs[i], block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                       genotypeDim=genoDim, verbose=FALSE)
      
      gds <- GdsGenotypeReader(gdsfile)
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      
      gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                            input.files=c(files[i], markfile, posfile), chromosome,
                            input.type="MaCH", 
                            input.dosage=inputs[i], block.size=b)
      
      geno.orig <- getGenotype(genoData)
      
      close(genoData)
      
      # change a snp
      gds <- openfn.gds(gdsfile, readonly=FALSE)
      orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
      val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
      write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
      closefn.gds(gds)
      
      # check exception
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      
      geno <- getGenotype(genoData)
      
      checkException(checkEquals(geno, geno.orig))
      checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                           input.files=c(files[i], markfile, posfile), chromosome,
                                           input.type="MaCH", 
                                           input.dosage=inputs[i], block.size=b))
      
      
      close(genoData)
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}


# tests adding only a subset of samples/snps
test_check_mach_subset <- function() {
  
  probfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlprob",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mldose",
                          package="GWASdata")
  markfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlinfo",
                          package="GWASdata")
  posfile <- system.file("extdata", "imputation", "MaCH", "mach1.snp.position",
                         package="GWASdata")
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  # set up scan.df for subsetting later.
  dosages <- read.table(dosefile, header=FALSE, stringsAsFactors=FALSE)
  # remove 5 random samples
  i_samp_rm <- sample(1:nrow(dosages), nrow(dosages)-5)
  scan.df <- data.frame(sampleID=c("MISSING->MISSING", dosages$V1[i_samp_rm]), stringsAsFactors=FALSE)
  scan.df$scanID <- 200:(200+nrow(scan.df)-1)
  
  markers <- read.table(markfile, header=TRUE, stringsAsFactors=FALSE)
  i_snp_rm <- sample(1:nrow(markers), 5)
  #message(paste(i_snp_rm, collapse=" "))
  #snp.names <- markers$SNP[i_snp_rm] # remove 5 random SNPs
  
  snp.id.start <- 100
  
  nsnp <- nrow(markers) - length(i_snp_rm)
  nscan <- nrow(scan.df)
  
  
  # random snp and scan to remove
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(which(scan.df$sampleID != "MISSING->MISSING"), 1)
  
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 1)
  genoDim <- "snp,scan"
  for (b in blocks) {
    for (i in 1:2) {
      
      # test reading snp.names and scan.df
      gdsImputedDosage(input.files=c(files[i], markfile, posfile), gds.filename=gdsfile, chromosome=22,
                       input.type="MaCH", input.dosage=inputs[i], block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                       scan.df=scan.df, snp.exclude=i_snp_rm, genotypeDim=genoDim,
                       snp.id.start=snp.id.start, verbose=FALSE)
      
      
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      
      gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                            input.files=c(files[i], markfile, posfile), chromosome,
                            input.type="MaCH", 
                            input.dosage=inputs[i], block.size=b,
                            snp.exclude=i_snp_rm)
      
      # change it in a random place
      close(genoData)
      
      gds <- openfn.gds(gdsfile, readonly=FALSE)
      orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp, i_scan), count=c(1,1))
      val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
      write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp, i_scan), count=c(1,1))
      closefn.gds(gds)
      
      # check exception
      
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                            input.files=c(files[i], markfile, posfile), chromosome,
                            input.type="MaCH", 
                            input.dosage=inputs[i], block.size=b,
                            snp.exclude=i_snp_rm))
      
      close(genoData)
      
      # change an added=FALSE sample
      tmp <- which(!scanAnnot$added)
      i_scan2 <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
      i_snp2 <- sample(1:nsnp, 1)
      gds <- openfn.gds(gdsfile, readonly=FALSE)
      # put the original value back
      write.gdsn(index.gdsn(gds, "genotype"), orig, start=c(i_snp,i_scan), count=c(1,1))
      # write 
      write.gdsn(index.gdsn(gds, "genotype"), 1, start=c(i_snp2, i_scan2), count=c(1,1))
      closefn.gds(gds)
      
      # check exception
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                           input.files=c(files[i], markfile, posfile), chromosome,
                                           input.type="MaCH", 
                                           input.dosage=inputs[i], block.size=b,
                                           snp.exclude=i_snp_rm))
      
      
      close(genoData)
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}




test_check_impute2 <- function() {
  probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                          package="GWASdata")
  sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                          package="GWASdata")
  
  samp <- read.table(sampfile, stringsAsFactors=FALSE, header=TRUE)
  samp <- samp[-1, ]
  
  dos <- read.table(probfile, header=FALSE, stringsAsFactors=FALSE)
  snps <- dos[, 2]
  
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  nsnp <- length(snps)
  nscan <- nrow(samp)
  
  # random snp and scan to remove
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(1:nscan, 1)

  # 33 lines in file
  blocks <- c(5000, 10, 32)
  genoDim <- "snp,scan"
  #for (genoDim in c("snp,scan", "scan,snp")) {
    for (b in blocks) {
        
      # make a normal one
      gdsImputedDosage(input.files=c(probfile, sampfile), gds.filename=gdsfile, chromosome=22,
                       input.type="IMPUTE2", input.dosage=FALSE, block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim,
                       verbose=FALSE)
      

      gds <- GdsGenotypeReader(gdsfile)
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      
      gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                            input.files=c(probfile, sampfile), chromosome,
                            input.type="IMPUTE2", 
                            input.dosage=FALSE, block.size=b)
      
      close(genoData)

      gds <- openfn.gds(gdsfile, readonly=FALSE)
      orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
      val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
      write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
      closefn.gds(gds)
      
      # check exception
      
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                           input.files=c(probfile, sampfile), chromosome,
                                           input.type="IMPUTE2", 
                                           input.dosage=FALSE, block.size=b))
      

      close(genoData)
      
    }
  
  unlink(c(gdsfile, snpfile, scanfile))
}


test_check_impute2_missing <- function() {
  probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                          package="GWASdata")
  sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                          package="GWASdata")
  
  samp <- read.table(sampfile, stringsAsFactors=FALSE, header=TRUE)
  samp <- samp[-1, ]
  
  dos <- read.table(probfile, header=FALSE, stringsAsFactors=FALSE)
  snps <- dos[, 2]
  
  # make a dosage missing
  newprobfile <- tempfile()
  dos[1, 6:8] <- -1
  write.table(dos, file=newprobfile, row.names=FALSE, col.names=FALSE)
  x <- read.table(newprobfile, stringsAsFactors=FALSE, header=FALSE)
  checkEquals(x, dos)
  
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  nsnp <- length(snps)
  nscan <- nrow(samp)
  
  # random snp and scan to remove
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(1:nscan, 1)
  
  # 33 lines in file
  blocks <- c(5000, 10, 32)
  genoDim <- "snp,scan"
  #for (genoDim in c("snp,scan", "scan,snp")) {
  for (b in blocks) {
    
    # make a normal one
    gdsImputedDosage(input.files=c(probfile, sampfile), gds.filename=gdsfile, chromosome=22,
                     input.type="IMPUTE2", input.dosage=FALSE, block.size=b,
                     snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim,
                     verbose=FALSE)
    
    
    gds <- GdsGenotypeReader(gdsfile)
    scanAnnot <- getobj(scanfile)
    snpAnnot <- getobj(snpfile)
    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    
    gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                          input.files=c(probfile, sampfile), chromosome,
                          input.type="IMPUTE2", 
                          input.dosage=FALSE, block.size=b)
    
    close(genoData)
    
    gds <- openfn.gds(gdsfile, readonly=FALSE)
    orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
    val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
    write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
    closefn.gds(gds)
    
    # check exception
    
    gds <- GdsGenotypeReader(gdsfile)
    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                         input.files=c(probfile, sampfile), chromosome,
                                         input.type="IMPUTE2", 
                                         input.dosage=FALSE, block.size=b))
    
    
    close(genoData)
    
  }
  
  unlink(c(gdsfile, snpfile, scanfile))
}



# tests adding only a subset of samples/snps
test_check_impute2_subset <- function() {
  probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                          package="GWASdata")
  sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                          package="GWASdata")
  
  samp <- read.table(sampfile, stringsAsFactors=FALSE, header=TRUE)
  samp <- samp[-1, ]
  i_samp_rm <- sample(-1:-nrow(samp), 5)
  scan.df <- data.frame(sampleID=paste(samp$ID_1[i_samp_rm], samp$ID_2[i_samp_rm]), stringsAsFactors=FALSE)
  scan.df <- rbind(data.frame(sampleID="MISSING MISSING", stringsAsFactors=FALSE), scan.df)
  scan.df$scanID <- 200:(200+nrow(scan.df)-1)
  
  dos <- read.table(probfile, header=FALSE, stringsAsFactors=FALSE)
  snps <- dos[, 2]
  # remove 5 random snps
  i_snp_rm <- sample(1:length(snps), 5)
  
  snp.id.start <- 100
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
  nsnp <- length(snps) - length(i_snp_rm)
  nscan <- nrow(samp)
  
  # random snp and scan to remove
  i_snp <- sample(1:nsnp, 1)
  i_scan <- sample(which(scan.df$sampleID != "MISSING MISSING"), 1)
  
  # 33 lines in file
  blocks <- c(5000, 10, 32, 1)
  genoDim <- "snp,scan"
  for (b in blocks) {
    
    # now the subset of samples/snps
    gdsImputedDosage(input.files=c(probfile, sampfile), gds.filename=gdsfile, chromosome=22,
                     input.type="IMPUTE2", input.dosage=FALSE, block.size=b,
                     snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                     scan.df=scan.df, snp.exclude=i_snp_rm, genotypeDim=genoDim,
                     snp.id.start=snp.id.start, verbose=FALSE)
    
    scanAnnot <- getobj(scanfile)
    snpAnnot <- getobj(snpfile)
    gds <- GdsGenotypeReader(gdsfile)
    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    
    gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                          input.files=c(probfile, sampfile), chromosome,
                          input.type="IMPUTE2", 
                          input.dosage=FALSE, block.size=b,
                          snp.exclude=i_snp_rm, verbose=FALSE)
        
    # change it in a random place
    close(genoData)
    
    gds <- openfn.gds(gdsfile, readonly=FALSE)
    orig <- read.gdsn(index.gdsn(gds, "genotype"), start=c(i_snp,i_scan), count=c(1,1))
    val <- ifelse(is.na(orig) | orig < 0, 1, orig-1)
    write.gdsn(index.gdsn(gds, "genotype"), val, start=c(i_snp,i_scan), count=c(1,1))
    closefn.gds(gds)
    
    # check exception
    
    gds <- GdsGenotypeReader(gdsfile)
    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                         input.files=c(probfile, sampfile), chromosome,
                                         input.type="IMPUTE2", 
                                         input.dosage=FALSE, block.size=b,
                                         snp.exclude=i_snp_rm, verbose=FALSE))
    
    close(genoData)
    
    # change an added=FALSE sample
    tmp <- which(!scanAnnot$added)
    i_scan2 <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
    i_snp2 <- sample(1:nsnp, 1)
    gds <- openfn.gds(gdsfile, readonly=FALSE)
    # put the original value back
    write.gdsn(index.gdsn(gds, "genotype"), orig, start=c(i_snp,i_scan), count=c(1,1))
    # write 
    write.gdsn(index.gdsn(gds, "genotype"), 1, start=c(i_snp2, i_scan2), count=c(1,1))
    closefn.gds(gds)
    
    # check exception
    gds <- GdsGenotypeReader(gdsfile)
    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    checkException(gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                                         input.files=c(probfile, sampfile), chromosome,
                                         input.type="IMPUTE2", 
                                         input.dosage=FALSE, block.size=b,
                                         snp.exclude=i_snp_rm))
    
    
    close(genoData)
  }
  unlink(c(gdsfile, snpfile, scanfile))
}
