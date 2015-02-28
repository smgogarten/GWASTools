library(logistf)

test_monomorphic <- function() {
    geno <- matrix(c(0,0,0,NA,
                     0,1,NA,2,
                     NA,2,2,2), nrow=3, byrow=TRUE)
    checkEquals(c(TRUE,FALSE,TRUE), GWASTools:::.monomorphic(geno, outcome="", model.type=""))

    ## case/control
    geno <- matrix(c(0,0,1,NA,
                     0,1,1,2,
                     NA,1,2,2), nrow=3, byrow=TRUE)
    checkEquals(c(TRUE,FALSE,TRUE), GWASTools:::.monomorphic(geno, outcome=c(0,0,1,1), model.type="logistic"))
}

.regModelData <- function(type=c("logistic", "linear", "poisson"), nsamp=100) {
    type <- match.arg(type)
    data.frame(outcome=switch(type,
                   logistic=rbinom(nsamp,1,0.4),
                   linear=rnorm(nsamp, mean=10, sd=2),
                   poisson=c(sample(1:20, 80, replace=TRUE), sample(21:40, 20, replace=TRUE))),
               covar=sample(letters[1:3], nsamp, replace=TRUE),
               genotype=sample(c(0,1,2), nsamp, replace=TRUE))
}

test_runReg_linear <- function() {
    mod <- "outcome ~ covar + genotype"
    dat <- .regModelData(type="linear")
    tmp <- GWASTools:::.runRegression(mod, dat, "linear",
                                      CI=0.95, robust=FALSE, LRtest=FALSE)
    checkIdentical(names(tmp), c("Est", "SE", "LL", "UL", "Wald.Stat", "Wald.pval"))
    checkTrue(all(!is.na(tmp)))

    fit <- lm(as.formula(mod), data=dat)
    checkEquals(summary(fit)$coef["genotype",1:2], tmp[1:2], checkNames=FALSE)
}

test_runReg_logistic <- function() {
    mod <- "outcome ~ covar + genotype"
    dat <- .regModelData(type="logistic")
    tmp <- GWASTools:::.runRegression(mod, dat, "logistic",
                                      CI=0.95, robust=FALSE, LRtest=FALSE)
    checkIdentical(names(tmp), c("Est", "SE", "LL", "UL", "Wald.Stat", "Wald.pval"))
    checkTrue(all(!is.na(tmp)))

    fit <- glm(as.formula(mod), data=dat, family=binomial())
    checkEquals(summary(fit)$coef["genotype",1:2], tmp[1:2], checkNames=FALSE)
}

test_runReg_poisson <- function() {
    mod <- "outcome ~ covar + genotype"
    dat <- .regModelData(type="poisson")
    tmp <- GWASTools:::.runRegression(mod, dat, "poisson",
                                      CI=0.95, robust=FALSE, LRtest=FALSE)
    checkIdentical(names(tmp), c("Est", "SE", "LL", "UL", "Wald.Stat", "Wald.pval"))
    checkTrue(all(!is.na(tmp)))

    fit <- glm(as.formula(mod), data=dat, family=poisson())
    checkEquals(summary(fit)$coef["genotype",1:2], tmp[1:2], checkNames=FALSE)
}

test_runReg_CI <- function() {
    mod <- "outcome ~ covar + genotype"
    dat <- .regModelData(type="linear")
    tmp <- GWASTools:::.runRegression(mod, dat, "linear",
                                      CI=0.95, robust=FALSE, LRtest=FALSE)
    tmp2 <- GWASTools:::.runRegression(mod, dat, "linear",
                                       CI=0.90, robust=FALSE, LRtest=FALSE)
    checkTrue(tmp["LL"] < tmp2["LL"])
    checkTrue(tmp["UL"] > tmp2["UL"])
}

test_runReg_robust <- function() {
    mod <- "outcome ~ covar + genotype"
    tmp <- GWASTools:::.runRegression(mod, .regModelData(type="linear"), "linear",
                                      CI=0.95, robust=TRUE, LRtest=FALSE)
    checkTrue(all(!is.na(tmp)))
}

test_runReg_LR <- function() {
    mod <- "outcome ~ covar + genotype"
    tmp <- GWASTools:::.runRegression(mod, .regModelData(type="linear"), "linear",
                                      CI=0.95, robust=FALSE, LRtest=TRUE)
    checkIdentical(names(tmp), c("Est", "SE", "LL", "UL", "Wald.Stat", "Wald.pval",
                                 "LR.Stat", "LR.pval"))
    checkTrue(all(!is.na(tmp)))
}

test_runReg_GxE <- function() {
    mod <- as.formula("outcome ~ covar + covar:genotype + genotype")
    tmp <- GWASTools:::.runRegression(mod, .regModelData(type="linear"), "linear",
                                      CI=0.95, robust=FALSE, LRtest=FALSE)
    checkIdentical(names(tmp), c("Est", "SE", "LL", "UL", "Wald.Stat", "Wald.pval",
                                 "GxE.Stat", "GxE.pval", "Joint.Stat", "Joint.pval"))
    checkTrue(all(!is.na(tmp)))
}
    
test_runFirth <- function() {
    mod <- "outcome ~ covar + genotype"
    dat <- .regModelData()
    tmp <- GWASTools:::.runFirth(mod, dat, CI=0.95, PPLtest=FALSE)
    checkIdentical(names(tmp), c("Est", "SE", "LL", "UL", "Wald.Stat", "Wald.pval"))
    fm <- logistf(as.formula(mod), data=dat, pl=FALSE)
    checkEquals(tmp["Est"], coef(fm)["genotype"], checkNames=FALSE)
    checkEquals(tmp["LL"], fm$ci.lower["genotype"], checkNames=FALSE)
    checkEquals(tmp["UL"], fm$ci.upper["genotype"], checkNames=FALSE)
    checkEquals(tmp["Wald.pval"], fm$prob["genotype"], checkNames=FALSE)
    checkEquals(GWASTools:::.waldTest(tmp["Est"], (tmp["SE"])^2),
                tmp[c("Wald.Stat", "Wald.pval")], checkNames=FALSE)
}

test_runFirth_CI <- function() {
    mod <- "outcome ~ covar + genotype"
    dat <- .regModelData(type="linear")
    tmp <- GWASTools:::.runFirth(mod, dat, CI=0.95, PPLtest=FALSE)
    checkEquals(GWASTools:::.CI(tmp["Est"], tmp["SE"], 0.95),
                tmp[c("LL", "UL")], checkNames=FALSE)
    
    tmp2 <- GWASTools:::.runFirth(mod, dat, CI=0.90, PPLtest=FALSE)
    checkTrue(tmp["LL"] < tmp2["LL"])
    checkTrue(tmp["UL"] > tmp2["UL"])
}

test_runFirth_PPL <- function() {
    mod <- "outcome ~ covar + genotype"
    dat <- .regModelData()
    ind <- which(colnames(model.matrix(as.formula(mod), dat)) == "genotype")
    tmp <- GWASTools:::.runFirth(mod, dat, CI=0.95, PPLtest=TRUE, geno.index=ind)
    checkIdentical(names(tmp), c("Est", "SE", "LL", "UL", "Wald.Stat", "Wald.pval",
                                 "PPL.Stat", "PPL.pval"))
    fm <- logistf(as.formula(mod), data=dat, pl=TRUE, plconf=ind)
    checkEquals(tmp["Est"], coef(fm)["genotype"], checkNames=FALSE)
    checkEquals(tmp["LL"], fm$ci.lower["genotype"], checkNames=FALSE)
    checkEquals(tmp["UL"], fm$ci.upper["genotype"], checkNames=FALSE)
    checkEquals(tmp["PPL.pval"], fm$prob["genotype"], checkNames=FALSE)
}


.regGenoData <- function(nsnp=100, nsamp=50, chromosome=NULL) {
    if (is.null(chromosome)) {
        chromosome <- rep(c(1L,23L), each=nsnp/2)
    } else {
        chromosome <- rep(as.integer(chromosome), nsnp)
    }
    geno <- matrix(sample(c(0,1,2,NA), nsnp*nsamp, replace=TRUE), nrow=nsnp, ncol=nsamp)
    geno[1,] <- 0 ## make one monomorphic
    mgr <- MatrixGenotypeReader(geno, snpID=1:nsnp, scanID=1:nsamp,
                                chromosome=chromosome,
                                position=1:nsnp)

    scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=1:nsamp,
      sex=sample(c("M","F"), nsamp, replace=TRUE),
      age=round(rnorm(nsamp, mean=40, sd=10)),
      site=sample(letters[1:3], nsamp, replace=TRUE),   
      status=rbinom(nsamp,1,0.4),
      trait=rnorm(nsamp, mean=10, sd=2)))

    GenotypeData(mgr, scanAnnot=scanAnnot)
}

test_blocks <- function() {
    outcome <- "trait"
    model.type <- "linear"

    genoData <- .regGenoData()
    assoc1 <- assocRegression(genoData,
                          outcome = outcome,
                          model.type = model.type)
    assoc2 <- assocRegression(genoData,
                          outcome = outcome,
                          model.type = model.type,
                          block.size=33)
    checkEquals(assoc1, assoc2)
}

test_snps <- function() {
    outcome <- "trait"
    model.type <- "linear"

    genoData <- .regGenoData()
    assoc1 <- assocRegression(genoData,
                          outcome = outcome,
                          model.type = model.type)
    assoc2a <- assocRegression(genoData,
                          outcome = outcome,
                          model.type = model.type,
                          snpStart=1, snpEnd=50)
    assoc2b <- assocRegression(genoData,
                          outcome = outcome,
                          model.type = model.type,
                          snpStart=51, snpEnd=100)
    checkEquals(assoc1, rbind(assoc2a, assoc2b))
}


.checkAssoc <- function(genoData, outcome, model.type, covar=NULL,
                        scan.exclude=NULL, robust=FALSE, LRtest=FALSE, ivar=NULL,
                        gene.action="additive") {

    if (!is.null(covar)) covar.list <- list(covar) else covar.list <- NULL  
    if (!is.null(ivar)) ivar.list <- list(ivar) else ivar.list <- NULL    
    assoc1 <- assocTestRegression(genoData,
                    outcome = outcome,
                    model.type = model.type,
                    covar.list = covar.list,
                    ivar.list = ivar.list,
                    gene.action.list = list(gene.action),
                    scan.exclude = scan.exclude,
                    robust = robust,
                    LRtest = LRtest)

    assoc2 <- assocRegression(genoData,
                    outcome = outcome,
                    model.type = model.type,
                    covar = covar,
                    ivar = ivar,
                    gene.action = gene.action,
                    scan.exclude = scan.exclude,
                    robust = robust,
                    LRtest = LRtest)

    cols2 <- intersect(c("snpID", "n", "MAF", "minor.allele", "Est", "SE",
                         unlist(lapply(c("Wald", "LR", "GxE", "Joint"), paste, c("Stat", "pval"), sep="."))),
                       names(assoc2))
    assoc2 <- assoc2[,cols2]
    cols1 <- c(cols2[1], paste0("model.1.", c(cols2[2:4], paste(gene.action, cols2[5:length(cols2)], "G", sep="."))))
    if (!is.null(ivar)) {
        for (x in c("Stat", "pval")) {
            cols1 <- sub(paste0("GxE.", x, ".G"), paste0("Wald.", x, ".G:", ivar), cols1)
            cols1 <- sub(paste0("Joint.", x, ".G"), paste0("Wald.", x, ".G.Joint"), cols1)
        }
    }
    assoc1 <- assoc1[,cols1]
    names(assoc1) <- cols2

    checkEquals(assoc1, assoc2)
}


test_logistic <- function() {
    outcome <- "status"
    model.type <- "logistic"
    covar <- c("age","sex")

    genoData <- .regGenoData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkAssoc(genoData, outcome, model.type, covar, scan.exclude)
}

test_linear <- function() {
    outcome <- "trait"
    model.type <- "linear"
    covar <- c("age","sex")

    genoData <- .regGenoData()
    scan.exclude <- sample(getScanID(genoData), 10)

    .checkAssoc(genoData, outcome, model.type, covar, scan.exclude)
}

test_robust <- function() {
    outcome <- "trait"
    model.type <- "linear"
    covar <- c("age","sex")

    genoData <- .regGenoData()

    .checkAssoc(genoData, outcome, model.type, covar, robust=TRUE)
}

test_LR <- function() {
    outcome <- "trait"
    model.type <- "linear"
    covar <- c("age","sex")

    genoData <- .regGenoData()

    .checkAssoc(genoData, outcome, model.type, covar, LRtest=TRUE)
}

test_GxE <- function() {
    outcome <- "trait"
    model.type <- "linear"
    covar <- c("age","sex")
    ivar <- "sex"

    genoData <- .regGenoData()

    .checkAssoc(genoData, outcome, model.type, covar, ivar=ivar)
}

test_gene.action <- function() {  
    outcome <- "trait"
    model.type <- "linear"
    covar <- c("age","sex")

    genoData <- .regGenoData()

    .checkAssoc(genoData, outcome, model.type, covar, gene.action="recessive")
    .checkAssoc(genoData, outcome, model.type, covar, gene.action="dominant")
}

## for now, just check there are no errors
test_firth <- function() {
    outcome <- "status"
    model.type <- "firth"
    covar <- c("age","sex","site")
    PPLtest <- TRUE

    genoData <- .regGenoData()
    scan.exclude <- sample(getScanID(genoData), 10)

    assoc <- assocRegression(genoData,
                          outcome = outcome,
                          model.type = model.type,
                          covar = covar,
                          scan.exclude = scan.exclude,
                          PPLtest = PPLtest)
    checkTrue(!all(is.na(assoc$Est)))
}

test_effectAllele <- function() {
    outcome <- "trait"
    model.type <- "linear"
    covar <- c("age","sex")

    genoData <- .regGenoData()
    assoc1 <- assocRegression(genoData,
                           outcome = outcome,
                           model.type = model.type,
                           covar = covar,
                           effectAllele="minor")
    assoc2 <- assocRegression(genoData,
                           outcome = outcome,
                           model.type = model.type,
                           covar = covar,
                           effectAllele="alleleA")

    Bmin <- assoc1$minor.allele == "B"
    checkEquals(assoc1$Est[Bmin], -assoc2$Est[Bmin])
    checkEquals(assoc1$Est[!Bmin], assoc2$Est[!Bmin])
}

test_Ychr <- function() {
    genoData <- .regGenoData(chromosome=25)
    scan.exclude <- getScanID(genoData)[getSex(genoData) == "F"]

    outcome <- "trait"
    model.type <- "linear"
    assoc1 <- assocRegression(genoData,
                    outcome = outcome,
                    model.type = model.type)
    assoc2 <- assocRegression(genoData,
                    outcome = outcome,
                    model.type = model.type,
                    scan.exclude=scan.exclude)
    checkEquals(assoc1, assoc2)

    checkException(assocRegression(genoData,
                    outcome = outcome, covar="sex",
                    model.type = model.type))

    .checkAssoc(genoData, outcome, model.type, scan.exclude=scan.exclude)
}
