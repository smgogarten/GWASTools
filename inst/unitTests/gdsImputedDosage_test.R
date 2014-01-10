test_probToDosage_beagle <- function() {
  probfile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.gprobs",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.dose",
                          package="GWASdata")
  
  prob <- read.table(probfile, as.is=TRUE, header=TRUE)
  prob <- as.matrix(prob[,4:ncol(prob)])
  
  dose <- read.table(dosefile, as.is=TRUE, header=TRUE)
  dose <- 2 - as.matrix(dose[,4:ncol(dose)])
  
  checkEquals(dose, GWASTools:::.probToDosage(prob, BB=TRUE), tolerance=0.0001)
}

test_probToDosage_mach <- function() {
  probfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlprob",
                          package="GWASdata")
  dosefile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mldose",
                          package="GWASdata")
  
  prob <- read.table(probfile, as.is=TRUE, header=FALSE)
  prob <- as.matrix(prob[,3:ncol(prob)])
  
  dose <- read.table(dosefile, as.is=TRUE, header=FALSE)
  dose <- as.matrix(dose[,3:ncol(dose)])
  
  test <- GWASTools:::.probToDosage(prob, BB=FALSE)
  # make colnames agree to pass check
  colnames(test) <- colnames(dose)
  checkEquals(dose, test, tolerance=0.001)
}

test_beagle <- function() {
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
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 40, 99)
  for (genoDim in c("snp,scan", "scan,snp")) {
    for (b in blocks) {
      for (i in 1:2) {
        gdsImputedDosage(input.files=c(files[i], markfile), gds.filename=gdsfile, chromosome=22,
                         input.type="BEAGLE", input.dosage=inputs[i], block.size=b,
                         snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim)
        
        gds <- GdsGenotypeReader(gdsfile)
        scanAnnot <- getobj(scanfile)
        snpAnnot <- getobj(snpfile)
        genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
        geno <- getGenotype(genoData)
        #alleleA <- getVariable(genoData, "alleleA")
        #alleleB <- getVariable(genoData, "alleleB")
        alleleA <- getAlleleA(genoData)
        alleleB <- getAlleleB(genoData)
        checkIdentical(snpAnnot$alleleA, alleleA)
        checkIdentical(snpAnnot$alleleB, alleleB)
        
        dat <- read.table(dosefile, as.is=TRUE, header=TRUE)
        dose <- 2 - as.matrix(dat[,4:ncol(dat)])
        dimnames(dose) <- NULL
        checkEquals(dose, geno, tolerance=0.0001)
        checkIdentical(names(dat)[-1:-3], scanAnnot$sampleID)
        
        mark <- read.table(markfile, as.is=TRUE, header=FALSE)
        checkIdentical(mark[,1], snpAnnot$snp)
        checkIdentical(mark[,2], snpAnnot$position)
        checkIdentical(mark[,3], snpAnnot$alleleA)
        checkIdentical(mark[,4], snpAnnot$alleleB)
        
        close(genoData)
      }
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}


test_beagle_missing <- function() {
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
  
  files <- c(newprobfile, newdosefile)
  inputs <- c(FALSE, TRUE)
  # 100 lines in file
  blocks <- c(5000, 40, 99)
  for (genoDim in c("snp,scan", "scan,snp")) {
    for (b in blocks) {
      for (i in 1:2) {
        gdsImputedDosage(input.files=c(files[i], markfile), gds.filename=gdsfile, chromosome=22,
                         input.type="BEAGLE", input.dosage=inputs[i], block.size=b,
                         snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim)
        
        gds <- GdsGenotypeReader(gdsfile)
        scanAnnot <- getobj(scanfile)
        snpAnnot <- getobj(snpfile)
        genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
        geno <- getGenotype(genoData)
        
        dat <- read.table(newdosefile, as.is=TRUE, header=TRUE)
        dose <- 2 - as.matrix(dat[,4:ncol(dat)])
        dose[dose < 0 | dose > 2] <- NA
        dimnames(dose) <- NULL
        checkEquals(dose, geno, tolerance=0.0001)
        
        close(genoData)
      }
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}


# tests adding only a subset of samples/snps
test_beagle_subset <- function() {
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
  
  snp.id.start <- 100
  
  gdsfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()
  
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
                       snp.id.start=snp.id.start)
      
      dat <- read.table(dosefile, as.is=TRUE, header=TRUE)
      dose <- 2 - as.matrix(dat[,4:ncol(dat)])
      dimnames(dose) <- NULL
      # samples were switched and a missing sample inserted in the middle
      # five random SNPs not included
      dose <- cbind(dose[-i_snp_rm, 2], rep(NA, nrow(dose)-length(i_snp_rm)), dose[-i_snp_rm, 1])
      gds <- GdsGenotypeReader(gdsfile)
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      geno <- getGenotype(genoData)
      alleleA <- getAlleleA(genoData)
      alleleB <- getAlleleB(genoData)
      checkIdentical(snpAnnot$alleleA, alleleA)
      checkIdentical(snpAnnot$alleleB, alleleB)
      
      checkEquals(dose, geno, tolerance=0.0001)
      #sampIDs <- c(names(dat)[5], "missing", names(dat)[4])
      checkIdentical(scan.df$sampleID, scanAnnot$sampleID)
      checkEquals(scanAnnot$scanID, scan.df$scanID)
      
      mark <- read.table(markfile, as.is=TRUE, header=FALSE)
      checkIdentical(mark[-i_snp_rm, 1], snpAnnot$snp)
      checkIdentical(mark[-i_snp_rm, 2], snpAnnot$position)
      checkIdentical(mark[-i_snp_rm, 3], snpAnnot$alleleA)
      checkIdentical(mark[-i_snp_rm, 4], snpAnnot$alleleB)
      
      checkIdentical(1:nsnp(genoData) + as.integer(snp.id.start-1), getSnpID(genoData))
      
      close(genoData)
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}


test_mach <- function() {
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
  
  dosages <- read.table(dosefile, header=FALSE, stringsAsFactors=FALSE)
  
  markers <- read.table(markfile, header=TRUE, stringsAsFactors=FALSE)
  
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 500 lines in file
  blocks <- c(5000, 200, 499)
  for (genoDim in c("snp,scan", "scan,snp")) {
    for (b in blocks) {
      for (i in 1:2) {
        gdsImputedDosage(input.files=c(files[i], markfile, posfile), gds.filename=gdsfile, chromosome=22,
                         input.type="MaCH", input.dosage=inputs[i], block.size=b,
                         snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                         genotypeDim=genoDim)
        
        gds <- GdsGenotypeReader(gdsfile)
        scanAnnot <- getobj(scanfile)
        snpAnnot <- getobj(snpfile)
        genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
        geno <- getGenotype(genoData)
        # alleleA <- getVariable(genoData, "alleleA")
        # alleleB <- getVariable(genoData, "alleleB")
        alleleA <- getAlleleA(genoData)
        alleleB <- getAlleleB(genoData)
        checkIdentical(snpAnnot$alleleA, alleleA)
        checkIdentical(snpAnnot$alleleB, alleleB)
        
        dat <- read.table(dosefile, as.is=TRUE, header=FALSE)
        #samples <- as.data.frame(matrix(unlist(strsplit(dat[,1], "->")), ncol=2, byrow=TRUE))
        checkIdentical(scanAnnot$sampleID, dat[,1])
        #checkIdentical(scanAnnot$ID_2, dat[,1])
        dose <-  t(as.matrix(dat[,3:ncol(dat)]))
        dimnames(dose) <- NULL
        checkEquals(dose, geno, tolerance=0.001)
        
        mark <- read.table(markfile, as.is=TRUE, header=TRUE)
        checkIdentical(mark[,1], snpAnnot$snp)
        checkIdentical(mark[,2], snpAnnot$alleleA)
        checkIdentical(mark[,3], snpAnnot$alleleB)
        
        pos <- read.table(posfile, as.is=TRUE, header=TRUE)
        checkIdentical(pos[,1], snpAnnot$snp)
        checkIdentical(pos[,2], snpAnnot$position)
        
        close(genoData)
      }
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}


test_mach_missing <- function() {
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
  
  
  files <- c(newprobfile, newdosefile)
  inputs <- c(FALSE, TRUE)
  # 500 lines in file
  blocks <- c(5000, 1)
  genoDim <- c("snp,scan")
  for (b in blocks) {
    for (i in 1:2) {
      gdsImputedDosage(input.files=c(files[i], markfile, posfile), gds.filename=gdsfile, chromosome=22,
                       input.type="MaCH", input.dosage=inputs[i], block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                       genotypeDim=genoDim)
      
      gds <- GdsGenotypeReader(gdsfile)
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      geno <- getGenotype(genoData)
      
      dat <- read.table(newdosefile, as.is=TRUE, header=FALSE)
      dose <-  t(as.matrix(dat[,3:ncol(dat)]))
      dose[dose < 0 | dose > 2] <- NA
      dimnames(dose) <- NULL
      checkEquals(dose, geno, tolerance=0.001)
      
      close(genoData)
      
    }
  }
  
  unlink(c(gdsfile, snpfile, scanfile, newprobfile, newdosefile))
}


# tests adding only a subset of samples/snps
test_mach_subset <- function() {
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
  
  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  # 500 lines in file
  blocks <- c(5000, 200, 499, 1)
  genoDim <- "snp,scan"
  for (b in blocks) {
    for (i in 1:2) {
      # test scan.df and snp.names
      gdsImputedDosage(input.files=c(files[i], markfile, posfile), gds.filename=gdsfile,
                       chromosome=22,
                       input.type="MaCH", input.dosage=inputs[i], block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                       genotypeDim=genoDim,
                       scan.df=scan.df, snp.exclude=i_snp_rm,
                       snp.id.start=snp.id.start)
      
      #dose <- cbind(dose[i_snp_rm, 2], rep(NA, nrow(dose)-length(i_snp_rm)), dose[i_snp_rm, 1])
      
      gds <- GdsGenotypeReader(gdsfile)
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      geno <- getGenotype(genoData)
      alleleA <- getAlleleA(genoData)
      alleleB <- getAlleleB(genoData)
      checkIdentical(snpAnnot$alleleA, alleleA)
      checkIdentical(snpAnnot$alleleB, alleleB)
      
      dat <- read.table(dosefile, as.is=TRUE, header=FALSE, stringsAsFactors=FALSE)
      checkIdentical(scanAnnot$sampleID, scan.df$sampleID)
      dose <-  t(as.matrix(dat[,3:ncol(dat)]))
      # 5 random samples were removed and a missing sample was inserted at the top
      # first SNP was not included
      dose <- dose[-i_snp_rm, i_samp_rm]
      dose <- cbind(rep(NA, nrow(dose)), dose)
      dimnames(dose) <- NULL
      checkEquals(dose, geno, tolerance=0.001)
      
      checkIdentical(scan.df$sampleID, scanAnnot$sampleID)
      checkEquals(scanAnnot$scanID, scan.df$scanID)
      
      mark <- read.table(markfile, as.is=TRUE, header=TRUE)
      checkIdentical(mark[-i_snp_rm,1], snpAnnot$snp)
      checkIdentical(mark[-i_snp_rm,2], snpAnnot$alleleA)
      checkIdentical(mark[-i_snp_rm,3], snpAnnot$alleleB)
      
      pos <- read.table(posfile, as.is=TRUE, header=TRUE)
      checkIdentical(pos[-i_snp_rm, 1], snpAnnot$snp)
      checkIdentical(pos[-i_snp_rm, 2], snpAnnot$position)

      checkIdentical(1:nsnp(genoData) + as.integer(snp.id.start-1), getSnpID(genoData))
      
      close(genoData)
      
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}



test_impute2 <- function() {
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
  
  # 33 lines in file
  blocks <- c(5000, 10, 32)
  for (genoDim in c("snp,scan", "scan,snp")) {
    for (b in blocks) {
      # make a normal one
      gdsImputedDosage(input.files=c(probfile, sampfile), gds.filename=gdsfile, chromosome=22,
                       input.type="IMPUTE2", input.dosage=FALSE, block.size=b,
                       snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim)
      
      gds <- GdsGenotypeReader(gdsfile)
      scanAnnot <- getobj(scanfile)
      snpAnnot <- getobj(snpfile)
      genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
      geno <- getGenotype(genoData)
      # alleleA <- getVariable(genoData, "alleleA")
      # alleleB <- getVariable(genoData, "alleleB")
      alleleA <- getAlleleA(genoData)
      alleleB <- getAlleleB(genoData)
      checkIdentical(snpAnnot$alleleA, alleleA)
      checkIdentical(snpAnnot$alleleB, alleleB)
      
      dat <- read.table(probfile, as.is=TRUE, header=FALSE)
      dose <- GWASTools:::.probToDosage(as.matrix(dat[,6:ncol(dat)]))
      dimnames(dose) <- NULL
      checkEquals(dose, geno, tolerance=0.0001)
      checkIdentical(dat[,1], snpAnnot$snp)
      checkIdentical(dat[,2], snpAnnot$rsID)
      checkIdentical(dat[,3], snpAnnot$position)
      checkIdentical(dat[,4], snpAnnot$alleleA)
      checkIdentical(dat[,5], snpAnnot$alleleB)
      
      samp <- read.table(sampfile, as.is=TRUE, header=FALSE, skip=2)
      checkIdentical(paste(samp[,1], samp[,2]), scanAnnot$sampleID)
      
      close(genoData)
      
    }
  }
  unlink(c(gdsfile, snpfile, scanfile))
}

test_impute2_missing <- function() {
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
  
  # 33 lines in file
  genoDim <- "snp,scan"
  blocks <- c(5000, 1)
  for (b in blocks) {
    # make a normal one
    gdsImputedDosage(input.files=c(newprobfile, sampfile), gds.filename=gdsfile, chromosome=22,
                     input.type="IMPUTE2", input.dosage=FALSE, block.size=b,
                     snp.annot.filename=snpfile, scan.annot.filename=scanfile, genotypeDim=genoDim)
    
    gds <- GdsGenotypeReader(gdsfile)
    scanAnnot <- getobj(scanfile)
    snpAnnot <- getobj(snpfile)
    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    geno <- getGenotype(genoData)
    
    dat <- read.table(newprobfile, as.is=TRUE, header=FALSE)
    dose <- GWASTools:::.probToDosage(as.matrix(dat[,6:ncol(dat)]))
    dose[dose < 0 | dose > 2] <- NA
    dimnames(dose) <- NULL
    checkEquals(dose, geno, tolerance=0.0001)
    
    close(genoData)
    
  }

  unlink(c(gdsfile, snpfile, scanfile, newprobfile))
}



# tests adding only a subset of samples/snps
test_impute2_subset <- function() {
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
  
  # 33 lines in file
  blocks <- c(5000, 10, 32, 1)
  genoDim <- "snp,scan"
  for (b in blocks) {
    
    # now the subset of samples/snps
    gdsImputedDosage(input.files=c(probfile, sampfile), gds.filename=gdsfile, chromosome=22,
                     input.type="IMPUTE2", input.dosage=FALSE, block.size=b,
                     snp.annot.filename=snpfile, scan.annot.filename=scanfile,
                     scan.df=scan.df, snp.exclude=i_snp_rm, genotypeDim=genoDim,
                     snp.id.start=snp.id.start)
    
    gds <- GdsGenotypeReader(gdsfile)
    scanAnnot <- getobj(scanfile)
    snpAnnot <- getobj(snpfile)
    genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    geno <- getGenotype(genoData)
    # alleleA <- getVariable(genoData, "alleleA")
    # alleleB <- getVariable(genoData, "alleleB")
    alleleA <- getAlleleA(genoData)
    alleleB <- getAlleleB(genoData)
    checkIdentical(snpAnnot$alleleA, alleleA)
    checkIdentical(snpAnnot$alleleB, alleleB)
    
    dat <- read.table(probfile, as.is=TRUE, header=FALSE)
    dose <- GWASTools:::.probToDosage(as.matrix(dat[,6:ncol(dat)]))
    dimnames(dose) <- NULL
    dose <- dose[-i_snp_rm, i_samp_rm]
    dose <- cbind(rep(NA, nrow(dose)), dose)
    checkEquals(dose, geno, tolerance=0.0001)
    checkIdentical(dat[-i_snp_rm, 1], snpAnnot$snp)
    checkIdentical(dat[-i_snp_rm, 2], snpAnnot$rsID)
    checkIdentical(dat[-i_snp_rm, 3], snpAnnot$position)
    checkIdentical(dat[-i_snp_rm, 4], snpAnnot$alleleA)
    checkIdentical(dat[-i_snp_rm, 5], snpAnnot$alleleB)
    
    samp <- read.table(sampfile, as.is=TRUE, header=FALSE, skip=2)
    #checkIdentical(paste(samp[, 1], samp[, 2]), scanAnnot$sampleID)
    checkIdentical(scan.df$sampleID, scanAnnot$sampleID)
    checkEquals(scan.df$scanID, scanAnnot$scanID)
    checkEquals(scan.df$scanID, getScanID(genoData))
    
    checkIdentical(1:nsnp(genoData) + as.integer(snp.id.start-1), getSnpID(genoData))
    
    close(genoData)
  }
  unlink(c(gdsfile, snpfile, scanfile))
}
