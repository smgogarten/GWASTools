library(survival)

.cphModelData <- function(nsamp=100) {
    set.seed(10); event <- rbinom(nsamp,1,0.4)
    set.seed(11); time.to.event <- rnorm(nsamp,mean=100,sd=10)
    set.seed(12); covar <- sample(letters[1:3], nsamp, replace=TRUE)
    set.seed(13); genotype <- sample(c(0,1,2), nsamp, replace=TRUE)
    data.frame(event, time.to.event, covar, genotype)
}

test_runCPH <- function() {
    dat <- .cphModelData()
    mod <- "surv ~ covar + genotype"
    tmp <- GWASTools:::.runCPH(mod, dat, "event", "time.to.event")
    checkIdentical(names(tmp), c("Est", "SE", "Wald.Stat", "Wald.pval"))
    checkTrue(all(!is.na(tmp)))

    surv <- Surv(dat$time.to.event, dat$event)
    cph <- coxph(as.formula(mod), data=dat)
    checkEquals(summary(cph)$coef["genotype",c(1,3,5)], tmp[c(1,2,4)], checkNames=FALSE)
}

test_runCPH_GxE <- function() {
    dat <- .cphModelData()
    mod <- "surv ~ covar + covar:genotype + genotype"
    tmp <- GWASTools:::.runCPH(mod, dat, "event", "time.to.event")
    checkTrue(all(c("GxE.Stat", "GxE.pval") %in% names(tmp)))
    checkTrue(!is.na(tmp["GxE.Stat"]))
}


.cphGenoData <- function(nsnp=100, nsamp=50, Xchrom=TRUE) {
    set.seed(14); geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE), nrow=nsnp, ncol=nsamp)
    geno[1,] <- 0 ## make one monomorphic
    if (Xchrom) chr <- rep(c(1L,23L), each=nsnp/2) else chr <- rep(1L, nsnp)
    mgr <- MatrixGenotypeReader(geno, snpID=1:nsnp, scanID=1:nsamp,
                                chromosome=chr, position=1:nsnp)

    set.seed(15); sex <- sample(c("M","F"), nsamp, replace=TRUE)
    set.seed(16); age <- round(rnorm(nsamp, mean=40, sd=10))
    set.seed(17); site <- sample(letters[1:3], nsamp, replace=TRUE)
    set.seed(18); event <- rbinom(nsamp,1,0.4)
    set.seed(19); time.to.event <- rnorm(nsamp,mean=100,sd=10)
    scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=1:nsamp,
      sex, age, site, event, time.to.event))

    GenotypeData(mgr, scanAnnot=scanAnnot)
}


test_CPH <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- c("age","sex","site")
    
    genoData <- .cphGenoData()
    set.seed(20); scan.exclude <- sample(getScanID(genoData), 10)

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

test_cluster <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- "age"
    cluster <- "site"
    
    genoData <- .cphGenoData()
    assoc <- assocCoxPH(genoData,
                       event,
                       time.to.event,
                       covar = covar,
                       cluster = cluster)
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
                       ivar = ivar,
                       LRtest = TRUE)
    checkTrue(!all(is.na(assoc$GxE.Stat)))
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
    
    checkEquals(with(assoc, 2*MAF*(1-MAF)*n.events > 75), assoc$maf.filter)
}

test_LR <- function() {
    event <- "event"
    time.to.event <- "time.to.event"
    covar <- c("age","sex","site")
    
    genoData <- .cphGenoData()

    assoc <- assocCoxPH(genoData,
                       event,
                       time.to.event,
                       covar = covar,
                       LRtest = TRUE)
    checkTrue(all(c("LR.Stat", "LR.pval") %in% names(assoc)))
    checkTrue(!all(is.na(assoc$LR.Stat)))
    # check that LR is generated for some tests where Wald had a warning
    checkTrue(!all(is.na(assoc$Wald.Stat) & is.na(assoc$LR.Stat)))
}
