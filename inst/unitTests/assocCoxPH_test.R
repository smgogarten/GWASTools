library(survival)

.cphModelData <- function(nsamp=100) {
    data.frame(event=rbinom(nsamp,1,0.4),
               time.to.event=rnorm(nsamp,mean=100,sd=10),
               covar=sample(letters[1:3], nsamp, replace=TRUE),
               genotype=sample(c(0,1,2), nsamp, replace=TRUE))
}

test_runCPH <- function() {
    dat <- .cphModelData()
    mod <- "surv ~ covar + genotype"
    tmp <- GWASTools:::.runCPH(mod, dat, "event", "time.to.event")
    checkIdentical(names(tmp), c("Est", "SE", "z.Stat", "z.pval"))
    checkTrue(all(!is.na(tmp)))

    surv <- Surv(dat$time.to.event, dat$event)
    cph <- coxph(as.formula(mod), data=dat)
    checkEquals(summary(cph)$coef["genotype",-2], tmp, checkNames=FALSE)
}

test_runCPH_GxE <- function() {
    dat <- .cphModelData()
    mod <- "surv ~ covar + covar:genotype + genotype"
    tmp <- GWASTools:::.runCPH(mod, dat, "event", "time.to.event")
    checkIdentical(names(tmp), c("Est", "SE", "z.Stat", "z.pval", "GxE.Stat", "GxE.pval"))
    checkTrue(all(!is.na(tmp)))
}


.cphGenoData <- function(nsnp=100, nsamp=50) {
    geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE), nrow=nsnp, ncol=nsamp)
    geno[1,] <- 0 ## make one monomorphic
    mgr <- MatrixGenotypeReader(geno, snpID=1:nsnp, scanID=1:nsamp,
                                chromosome=rep(c(1L,23L), each=nsnp/2),
                                position=1:nsnp)

    scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=1:nsamp,
      sex=sample(c("M","F"), nsamp, replace=TRUE),
      age=round(rnorm(nsamp, mean=40, sd=10)),
      site=sample(letters[1:3], nsamp, replace=TRUE),
      event=rbinom(nsamp,1,0.4),
      time.to.event=rnorm(nsamp,mean=100,sd=10)))

    GenotypeData(mgr, scanAnnot=scanAnnot)
}


.checkCPH <- function(genoData, event, time.to.event, covar,
                      scan.exclude=NULL, ivar=NULL, strata=NULL) {

    assoc1 <- assocTestCPH(genoData,
                       event,
                       time.to.event,
                       covars = covar,
                       GxE = ivar,
                       strata = strata,
                       scan.exclude = scan.exclude)
    
    assoc2 <- assocCoxPH(genoData,
                       event,
                       time.to.event,
                       covar = covar,
                       ivar = ivar,
                       strata = strata,
                       scan.exclude = scan.exclude,
                       effectAllele = "alleleA")

    if (is.null(ivar)) {
        cols2 <- intersect(c("snpID", "chr", "n.events", "MAF", "Est", "SE", "z.Stat", "z.pval"), names(assoc2))
        assoc2 <- assoc2[,cols2]

        assoc1$MAF <- ifelse(assoc1$chr == 23, assoc1$mafx, assoc1$maf)
        cols1 <- c(cols2[1:4], c("beta", "se", "z", "pval"))
        assoc1 <- assoc1[,cols1]
        names(assoc1) <- cols2
        
        checkEquals(assoc1, assoc2)
    } else {
        ## assoc1[[1]] is ignoring interaction term
        checkEquals(assoc1[[2]][,c("ge.lrtest", "ge.pval")], assoc2[,c("GxE.Stat", "GxE.pval")],
                    checkNames=FALSE)
    }
}

test_CPH <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- c("age","sex","site")
    
    genoData <- .cphGenoData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkCPH(genoData, event, time.to.event, covar, scan.exclude)
}

test_strata <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- "age"
    strata <- c("sex","site")
    
    genoData <- .cphGenoData()
    .checkCPH(genoData, event, time.to.event, covar, strata=strata)
}

test_GxE <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- c("age","sex","site")
    ivar <- "site"
    
    genoData <- .cphGenoData()
    .checkCPH(genoData, event, time.to.event, covar, ivar=ivar)
}
