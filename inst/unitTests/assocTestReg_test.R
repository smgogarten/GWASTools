.regSampleData <- function(nsnp=100, nsamp=50) {
    geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE), nrow=nsnp, ncol=nsamp)
    mgr <- MatrixGenotypeReader(geno, snpID=1:nsnp, scanID=1:nsamp,
                                chromosome=rep(c(1L,23L), each=nsnp/2),
                                position=1:nsnp)

    scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=1:nsamp,
      sex=sample(c("M","F"), nsamp, replace=TRUE),
      age=rnorm(nsamp, mean=40, sd=10),
      case.cntl.status=rbinom(nsamp,1,0.4)))
    scanAnnot$blood.pressure[scanAnnot$case.cntl.status==1] <- rnorm(sum(scanAnnot$case.cntl.status==1),mean=100,sd=10)
    scanAnnot$blood.pressure[scanAnnot$case.cntl.status==0] <- rnorm(sum(scanAnnot$case.cntl.status==0),mean=90,sd=5)

    GenotypeData(mgr, scanAnnot=scanAnnot)
}

test_blocks <- function() {
    outcome <- "blood.pressure"
    model.type <- "linear"

    genoData <- .regSampleData()
    assoc1 <- assocTestReg(genoData,
                          outcome = outcome,
                          model.type = model.type)
    assoc2 <- assocTestReg(genoData,
                          outcome = outcome,
                          model.type = model.type,
                          block.size=33)
    checkEquals(assoc1, assoc2)
}

test_snps <- function() {
    outcome <- "blood.pressure"
    model.type <- "linear"

    genoData <- .regSampleData()
    assoc1 <- assocTestReg(genoData,
                          outcome = outcome,
                          model.type = model.type)
    assoc2a <- assocTestReg(genoData,
                          outcome = outcome,
                          model.type = model.type,
                          snpStart=1, snpEnd=50)
    assoc2b <- assocTestReg(genoData,
                          outcome = outcome,
                          model.type = model.type,
                          snpStart=51, snpEnd=100)
    checkEquals(assoc1, rbind(assoc2a, assoc2b))
}


.checkAssoc <- function(genoData, outcome, model.type, covar.vec,
                        scan.exclude, robust, LRtest) {
    
    assoc1 <- assocTestRegression(genoData,
                    outcome = outcome,
                    model.type = model.type,
                    covar.list = list(covar.vec),
                    gene.action.list = "additive",
                    scan.exclude = scan.exclude,
                    robust = robust,
                    LRtest = LRtest)

    assoc2 <- assocTestReg(genoData,
                    outcome = outcome,
                    model.type = model.type,
                    covar.vec = covar.vec,
                    scan.exclude = scan.exclude,
                    robust = robust,
                    LRtest = LRtest)

    cols2 <- intersect(c("snpID", "n", "MAF", "minor.allele", "Est", "SE", "Wald.Stat", "Wald.pval", "LR.Stat", "LR.pval"), names(assoc2))
    assoc2 <- assoc2[,cols2]
    cols1 <- c(cols2[1], paste0("model.1.", c(cols2[2:4], paste0("additive.", cols2[5:length(cols2)], ".G"))))
    assoc1 <- assoc1[,cols1]
    names(assoc1) <- cols2

    ## for assocTestRegression, effect allele is minor.
    ## for assocTestReg, effect allele is A.
    Bmin <- assoc1$minor.allele %in% "B"
    assoc1$Est[Bmin] <- -1*assoc1$Est[Bmin]

    ## for MAF==0.5, definition of minor allele is opposite
    ## want to match assocTestMixedModel
    assoc1$minor.allele[assoc1$MAF %in% 0.5] <- "B"

    for (i in names(assoc2)) checkEquals(assoc1[,i], assoc2[,i])
}


test_logistic <- function() {
    outcome <- "case.cntl.status"
    model.type <- "logistic"
    covar.vec <- c("age","sex")
    robust <- FALSE
    LRtest <- FALSE

    genoData <- .regSampleData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkAssoc(genoData, outcome, model.type, covar.vec,
                scan.exclude, robust, LRtest)
}

test_linear <- function() {
    outcome <- "blood.pressure"
    model.type <- "linear"
    covar.vec <- c("age","sex")
    robust <- FALSE
    LRtest <- FALSE   

    genoData <- .regSampleData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkAssoc(genoData, outcome, model.type, covar.vec,
                scan.exclude, robust, LRtest)
}

test_robust <- function() {
    outcome <- "blood.pressure"
    model.type <- "linear"
    covar.vec <- c("age","sex")
    robust <- TRUE
    LRtest <- FALSE

    genoData <- .regSampleData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkAssoc(genoData, outcome, model.type, covar.vec,
                scan.exclude, robust, LRtest)
}

test_LR <- function() {
    outcome <- "blood.pressure"
    model.type <- "linear"
    covar.vec <- c("age","sex")
    robust <- FALSE
    LRtest <- TRUE

    genoData <- .regSampleData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkAssoc(genoData, outcome, model.type, covar.vec,
                scan.exclude, robust, LRtest)
}
