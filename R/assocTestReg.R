## get data frame with outcome and covariates
.modelData <- function(genoData, outcome, covar.vec) {
    mod.vars <- outcome
    if (!is.null(covar.vec)) {
        cvnames <- unique(unlist(strsplit(covar.vec,"[*:]")))
        mod.vars <- c(mod.vars, cvnames)
    }
    dat <- as.data.frame(getScanVariable(genoData, mod.vars))
    names(dat)[1] <- outcome ## fix name in case of only one column
    dat
}

## geno is a matrix of genotypes [sample,snp]
## chrom is chromosome corresponding to geno[,snp]
## keep is the sample index corresponding to geno[sample,]
.freqFromGeno <- function(genoData, geno, chrom, keep) {
    freq <- 0.5*colMeans(geno, na.rm=TRUE)
    ## for X chr
    if (XchromCode(genoData) %in% chrom) {
        ## which are on X chr
        Xidx <- chrom == XchromCode(genoData)
        ## males
        m <- (getSex(genoData) == "M")[keep]
        f <- (getSex(genoData) == "F")[keep]
        ## calculate allele freq for X
        freq[Xidx] <- (0.5 * colSums(geno[m,Xidx,drop=FALSE], na.rm=TRUE) +
                       colSums(geno[f,Xidx,drop=FALSE], na.rm=TRUE)) /
                           (colSums(!is.na(geno[m,Xidx,drop=FALSE])) +
                            2*colSums(!is.na(geno[f,Xidx,drop=FALSE])))
    }
    freq
}

## geno is a matrix of genotypes [sample,snp]
.monomorphic <- function(maf, geno, outcome, model.type) {
    mono <- is.na(maf) | maf == 0
    if (model.type %in% c("logistic", "firth")) {
        cc0 <- outcome == 0
        cc1 <- outcome == 1
        freq.0 <- 0.5*colMeans(geno[cc0,,drop=FALSE], na.rm=TRUE)
        freq.1 <- 0.5*colMeans(geno[cc1,,drop=FALSE], na.rm=TRUE)
        mono.cc0 <- is.na(freq.0) | freq.0 == 0 | freq.0 == 1
        mono.cc1 <- is.na(freq.1) | freq.1 == 0 | freq.1 == 1
        mono <- mono | mono.cc0 | mono.cc1
    }
    mono
}

.CI <- function(Est, SE, CI) {
    LL <- Est + qnorm((1-CI)*0.5)*SE
    UL <- Est + qnorm(1-((1-CI)*0.5))*SE
    c(LL=LL, UL=UL)
}

.waldTest <- function(Est, SE) {
    W <- (Est/SE)^2
    pval <- pchisq(W, df=1, lower.tail=FALSE)
    c(Wald.Stat=W, Wald.pval=pval)
}

## give the data frame a crazy name so we don't overwrite something important
.runRegression <- function(model.formula, ..model..data.., model.type, CI, robust, LRtest) {
    tryCatch({
        if (model.type == "linear") {
            mod <- lm(model.formula, data=..model..data..)
        } else if (model.type == "logistic") {
            mod <- glm(model.formula, data=..model..data.., family=binomial())
        }
        
        if (!robust) {
            ## model based variance
            Vhat <- vcov(mod)
        } else {
            ## sandwich variance
            Vhat <- vcovHC(mod, type="HC0")
        }
        Est <- unname(coef(mod)["genotype"])
        SE <- sqrt(Vhat["genotype","genotype"])
        ret <- c(Est=Est, SE=SE, .CI(Est, SE, CI), .waldTest(Est, SE))

        if (LRtest) {
            ## promote the data frame to the global environment so lrtest can find it
            ..model..data.. <<- ..model..data..
            mod2 <- lrtest(mod, "genotype")[2,]
            ret <- c(ret, LR.Stat=mod2[["Chisq"]], LR.pval=mod2[["Pr(>Chisq)"]])
        }
        ret
    }, warning=function(w) NA, error=function(e) NA)
}

.runFirth <- function(model.formula, model.data, CI, PPLtest) {
    tryCatch({
        ind <- which(names(model.data) == "genotype")
        mod <- logistf(model.formula, data=model.data, alpha=(1-CI),
                       pl=PPLtest, plconf=ind, dataout=FALSE)
        Est <- unname(coef(mod)[ind])
        SE <- sqrt(vcov(mod)[ind,ind])
        LL <- unname(mod$ci.lower[ind])
        UL <- unname(mod$ci.upper[ind])
        pval <- unname(mod$prob[ind])
        Stat <- qchisq(pval, df=1, lower.tail=FALSE)
        
        ret <- c(Est=Est, SE=SE, LL=LL, UL=UL)
        if (PPLtest) {
            c(ret, .waldTest(Est, SE), PPL.Stat=Stat, PPL.pval=pval)
        } else {
            c(ret, Wald.Stat=Stat, Wald.pval=pval)
        }
    }, warning=function(w) NA, error=function(e) NA)
}

    
assocTestReg <- function(genoData,
                         outcome,
                         model.type = c("linear", "logistic", "firth"),
                         covar.vec = NULL,
                         scan.exclude = NULL,
                         CI = 0.95,
                         robust = FALSE,
                         LRtest = FALSE,
                         PPLtest = FALSE,
                         snpStart = NULL,
                         snpEnd = NULL,
                         block.size = 5000,
                         verbose = TRUE) {

    ## check that test is valid
    model.type <- match.arg(model.type)
    
    ## set snpStart and snpEnd
    if (is.null(snpStart)) snpStart <- 1
    if (is.null(snpEnd)) snpEnd <- nsnp(genoData)

    ## set which samples to keep
    if (!is.null(scan.exclude)) {
        keep <- !(getScanID(genoData) %in% scan.exclude)
    } else {
        keep <- rep(TRUE, nscan(genoData))
    }

    ## get chromosome information
    chr <- getChromosome(genoData, index=snpStart:snpEnd)

    ## X chromosome check for sex variable
    if (XchromCode(genoData) %in% chr & !hasSex(genoData)) {
        stop("Sex values for the samples are required to compute MAF for X chromosome SNPs")
    }

    ## Y chromosome 
    if (YchromCode(genoData) %in% chr) {
        ## check for sex variable
        if (!hasSex(genoData)) {
            stop("Sex values for the samples are required for Y chromosome SNPs")
        }
        if (!all(chr == YchromCode(genoData))) {
            stop("Y chromosome must be analyzed separately")
        }
        ## only keep males
        keep <- keep & (getSex(genoData) == "M")
    }


    ## read in outcome and covariate data
    if (verbose) message("Reading in Phenotype and Covariate Data...")
    dat <- .modelData(genoData, outcome, covar.vec)
    ## identify samples with any missing data
    keep <- keep & complete.cases(dat)
    dat <- dat[keep,]
    model.formula <- as.formula(paste(paste(outcome,"~"), paste(c(covar.vec, "genotype"),collapse="+")))


    ## sample size, assuming no missing genotypes
    n <- nrow(dat)
    if (verbose) message("Running analysis with ", n, " Samples")

    ## number of SNPs in the segment
    nsnp.seg <- snpEnd - snpStart + 1
    ## determine number of SNP blocks
    nblocks <- ceiling(nsnp.seg/block.size)

    ## set up results matrix
    nv <- c("snpID","chr","n","MAF","minor.allele","Est","SE","LL","UL","Wald.Stat","Wald.pval")
    if (LRtest & model.type %in% c("linear", "logistic")) nv <- c(nv, "LR.Stat", "LR.pval")
    if (PPLtest & model.type == "firth") nv <- c(nv, "PPL.Stat", "PPL.pval")
    res <- matrix(NA, nrow=nsnp.seg, ncol=length(nv), dimnames=list(NULL, nv))
    reg.cols <- which(colnames(res) == "Est"):ncol(res)

    ## chromosome
    res[,"chr"] <- chr

    if (verbose) message("Beginning Calculations...")
    ## loop through blocks
    for (b in 1:nblocks) {
        
        ## keep track of time for rate reporting
        startTime <- Sys.time()
        
        snp.start.pos <- snpStart + (b-1)*block.size
        nsnp.block <- block.size
        if (snp.start.pos + nsnp.block > snpEnd) {
            nsnp.block <- snpEnd - snp.start.pos + 1
        }
        snp.end.pos <- snp.start.pos + nsnp.block - 1
        
        bidx <- ((b-1)*block.size+1):((b-1)*block.size+nsnp.block)
        
        ## get genotypes for the block
        geno <- getGenotype(genoData, snp=c(snp.start.pos, nsnp.block), scan=c(1,-1), drop=FALSE, transpose=TRUE)
        ## subset
        geno <- geno[keep, , drop=FALSE]
        
        ## allele frequency
        freq <- .freqFromGeno(genoData, geno, chr[bidx], keep)
        
        ## MAF
        maf <- ifelse(freq < 0.5, freq, 1-freq)
        res[bidx,"MAF"] <- maf
        ## minor allele coding:  A = 1, B = 0
        res[bidx,"minor.allele"] <- ifelse(freq < 0.5, 1, 0)

        ## check for monomorphic SNPs
        mono <- .monomorphic(maf, geno, dat[[outcome]], model.type)
        
        ## sample size
        n <- colSums(!is.na(geno))
        res[bidx, "n"] <- n

        ## loop through SNPs in block
        midx <- (1:nsnp.block)[!mono]
        for (i in midx) {
            mdat <- cbind(dat, genotype=geno[,i])
            mdat <- mdat[complete.cases(mdat),]
            if (model.type %in% c("linear", "logistic")) {
                tmp <- .runRegression(model.formula, mdat, model.type, CI, robust, LRtest)
            } else if (model.type == "firth") {
                tmp <- .runFirth(model.formula, mdat, CI, PPLtest)
            }
            res[bidx[i], reg.cols] <- tmp
        }
        
        endTime <- Sys.time()
        rate <- format(endTime - startTime, digits=4)
        
        if (verbose) message(paste("Block", b, "of", nblocks, "Completed -", rate))
    } ## end block loop

    ## results data frame
    res <- as.data.frame(res)

    ## add in snpID
    res$snpID <- getSnpID(genoData, index=snpStart:snpEnd)

    ## convert minor.allele coding back to A/B
    res[,"minor.allele"][res[,"minor.allele"] == 1] <- "A"
    res[,"minor.allele"][res[,"minor.allele"] == 0] <- "B"

    return(res)
}
