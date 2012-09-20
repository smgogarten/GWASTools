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

  ncfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()

  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  for (i in 1:2) {
    ncdfImputedDosage(input.files=c(files[i], markfile), ncdf.filename=ncfile, chromosome=22,
                      input.type="BEAGLE", input.dosage=inputs[i],
                      snp.annot.filename=snpfile, scan.annot.filename=scanfile)

    nc <- NcdfGenotypeReader(ncfile)
    scanAnnot <- getobj(scanfile)
    snpAnnot <- getobj(snpfile)
    genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    geno <- getGenotype(genoData)
    alleleA <- getVariable(genoData, "alleleA")
    alleleB <- getVariable(genoData, "alleleB")
    checkIdentical(snpAnnot$alleleA, alleleA)
    checkIdentical(snpAnnot$alleleB, alleleB)
  
    dat <- read.table(dosefile, as.is=TRUE, header=TRUE)
    dose <- 2 - as.matrix(dat[,4:ncol(dat)])
    dimnames(dose) <- NULL
    checkEquals(dose, geno, tolerance=0.0001)
    checkIdentical(names(dat)[-1:-3], scanAnnot$ID)

    mark <- read.table(markfile, as.is=TRUE, header=FALSE)
    checkIdentical(mark[,1], snpAnnot$marker)
    checkIdentical(mark[,2], snpAnnot$position)
    checkIdentical(mark[,3], snpAnnot$alleleA)
    checkIdentical(mark[,4], snpAnnot$alleleB)

    close(genoData)
  }
  
  unlink(c(ncfile, snpfile, scanfile))
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

  ncfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()

  files <- c(probfile, dosefile)
  inputs <- c(FALSE, TRUE)
  for (i in 1:2) {
    ncdfImputedDosage(input.files=c(files[i], markfile, posfile), ncdf.filename=ncfile, chromosome=22,
                      input.type="MaCH", input.dosage=inputs[i],
                      snp.annot.filename=snpfile, scan.annot.filename=scanfile)

    nc <- NcdfGenotypeReader(ncfile)
    scanAnnot <- getobj(scanfile)
    snpAnnot <- getobj(snpfile)
    genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
    geno <- getGenotype(genoData)
    alleleA <- getVariable(genoData, "alleleA")
    alleleB <- getVariable(genoData, "alleleB")
    checkIdentical(as.character(snpAnnot$alleleA), alleleA)
    checkIdentical(as.character(snpAnnot$alleleB), alleleB)
  
    dat <- read.table(dosefile, as.is=TRUE, header=FALSE)
    samples <- as.data.frame(matrix(unlist(strsplit(dat[,1], "->")), ncol=2, byrow=TRUE))
    checkIdentical(scanAnnot$ID_1, samples[,1])
    checkIdentical(scanAnnot$ID_2, samples[,2])
    dose <-  t(as.matrix(dat[,3:ncol(dat)]))
    dimnames(dose) <- NULL
    checkEquals(dose, geno, tolerance=0.001)

    mark <- read.table(markfile, as.is=TRUE, header=TRUE)
    checkIdentical(mark[,1], snpAnnot$SNP)
    checkIdentical(mark[,2], snpAnnot$alleleA)
    checkIdentical(mark[,3], snpAnnot$alleleB)

    pos <- read.table(posfile, as.is=TRUE, header=TRUE)
    checkIdentical(pos[,1], snpAnnot$SNP)
    checkIdentical(pos[,2], snpAnnot$position)

    close(genoData)
  }
  
  unlink(c(ncfile, snpfile, scanfile))
}

test_impute2 <- function() {
  probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                          package="GWASdata")
  sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                          package="GWASdata")

  ncfile <- tempfile()
  snpfile <- tempfile()
  scanfile <- tempfile()

  ncdfImputedDosage(input.files=c(probfile, sampfile), ncdf.filename=ncfile, chromosome=22,
                    input.type="IMPUTE2", input.dosage=FALSE,
                    snp.annot.filename=snpfile, scan.annot.filename=scanfile)

  nc <- NcdfGenotypeReader(ncfile)
  scanAnnot <- getobj(scanfile)
  snpAnnot <- getobj(snpfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  geno <- getGenotype(genoData)
  alleleA <- getVariable(genoData, "alleleA")
  alleleB <- getVariable(genoData, "alleleB")
  checkIdentical(snpAnnot$alleleA, alleleA)
  checkIdentical(snpAnnot$alleleB, alleleB)
  
  dat <- read.table(probfile, as.is=TRUE, header=FALSE)
  dose <- GWASTools:::.probToDosage(as.matrix(dat[,6:ncol(dat)]))
  dimnames(dose) <- NULL
  checkEquals(dose, geno, tolerance=0.0001)
  checkIdentical(dat[,1], snpAnnot$SNP)
  checkIdentical(dat[,2], snpAnnot$rsID)
  checkIdentical(dat[,3], snpAnnot$position)
  checkIdentical(dat[,4], snpAnnot$alleleA)
  checkIdentical(dat[,5], snpAnnot$alleleB)

  samp <- read.table(sampfile, as.is=TRUE, header=FALSE, skip=2)
  checkIdentical(samp[,1], scanAnnot$ID_1)
  checkIdentical(samp[,2], scanAnnot$ID_2)

  close(genoData)
  
  unlink(c(ncfile, snpfile, scanfile))
}
