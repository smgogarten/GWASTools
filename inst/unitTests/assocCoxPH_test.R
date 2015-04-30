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


.cphGenoData <- function(nsnp=100, nsamp=50, Xchrom=TRUE) {
    geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE), nrow=nsnp, ncol=nsamp)
    geno[1,] <- 0 ## make one monomorphic
    if (Xchrom) chr <- rep(c(1L,23L), each=nsnp/2) else chr <- rep(1L, nsnp)
    mgr <- MatrixGenotypeReader(geno, snpID=1:nsnp, scanID=1:nsamp,
                                chromosome=chr, position=1:nsnp)

    scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=1:nsamp,
      sex=sample(c("M","F"), nsamp, replace=TRUE),
      age=round(rnorm(nsamp, mean=40, sd=10)),
      site=sample(letters[1:3], nsamp, replace=TRUE),
      event=rbinom(nsamp,1,0.4),
      time.to.event=rnorm(nsamp,mean=100,sd=10)))

    GenotypeData(mgr, scanAnnot=scanAnnot)
}


test_CPH <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- c("age","sex","site")
    
    genoData <- .cphGenoData()
    scan.exclude <- sample(getScanID(genoData), 10)

    assoc <- assocCoxPH(genoData,
                       event,
                       time.to.event,
                       covar = covar,
                       scan.exclude = scan.exclude)
    checkTrue(!all(is.na(assoc$Est)))
}

test_strata <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- "age"
    strata <- c("sex","site")
    
    genoData <- .cphGenoData()
    assoc <- assocCoxPH(genoData,
                       event,
                       time.to.event,
                       covar = covar,
                       strata = strata)
    checkTrue(!all(is.na(assoc$Est)))
}

test_GxE <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- c("age","sex","site")
    ivar <- "site"
    
    genoData <- .cphGenoData()
    assoc <- assocCoxPH(genoData,
                       event,
                       time.to.event,
                       covar = covar,
                       ivar = ivar)
    checkTrue(!all(is.na(assoc$Est)))
}

test_filter <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- c("age","sex","site")
    
    genoData <- .cphGenoData(nsamp=500, Xchrom=FALSE)

    assoc <- assocCoxPH(genoData,
                       event,
                       time.to.event,
                       covar = covar)
    
    checkEquals(with(assoc, 2*MAF*(1-MAF)*n.events > 75), assoc$filter)
}
