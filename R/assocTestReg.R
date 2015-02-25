## get data frame with outcome and covariates
.modelData <- function(genoData, outcome, covar, ivar=NULL) {
    mod.vars <- outcome
    if (!is.null(covar)) {
        cvnames <- unique(unlist(strsplit(covar,"[*:]")))
        mod.vars <- c(mod.vars, cvnames)
        if (!all(cvnames %in% getScanVariableNames(genoData))) {
            stop("Not all variables in covar found in scan annotation of genoData")
        }
    }
    if (!is.null(ivar)) {
        if (is.null(covar) | !(ivar %in% cvnames)) {
            stop("ivar should also be present in covar")
        }
    }
    dat <- as.data.frame(getScanVariable(genoData, mod.vars))
    names(dat)[1] <- outcome ## fix name in case of only one column
    dat
}

## geno is a matrix of genotypes [snp,sample]
## chrom is chromosome corresponding to snp
## keep is the sample index corresponding to sample
.freqFromGeno <- function(genoData, geno, chrom, keep) {
    freq <- 0.5*rowMeans(geno, na.rm=TRUE)
    ## for X chr
    if (XchromCode(genoData) %in% chrom) {
        ## which are on X chr
        Xidx <- chrom == XchromCode(genoData)
        ## males
        m <- (getSex(genoData) == "M")[keep]
        f <- (getSex(genoData) == "F")[keep]
        ## calculate allele freq for X
        freq[Xidx] <- (0.5 * rowSums(geno[Xidx,m,drop=FALSE], na.rm=TRUE) +
                       rowSums(geno[Xidx,f,drop=FALSE], na.rm=TRUE)) /
                           (rowSums(!is.na(geno[Xidx,m,drop=FALSE])) +
                            2*rowSums(!is.na(geno[Xidx,f,drop=FALSE])))
    }
    freq
}

# A function to Transform Genotypes based on gene action
.transformGenotype <- function(geno, gene.action) {
    switch(gene.action,
           additive = {},
           dominant = {geno[geno == 1] <- 2},
           recessive = {geno[geno == 1] <- 0}
           )
    geno
}

## geno is a matrix of genotypes [snp,sample]
.monomorphic <- function(geno, outcome, model.type) {
    freq <- 0.5*rowMeans(geno, na.rm=TRUE)
    mono <- is.na(freq) | freq == 0 | freq == 1
    if (model.type %in% c("logistic", "firth")) {
        for (case.status in c(0,1)) {
            cc <- outcome == case.status
            freq <- 0.5*rowMeans(geno[,cc,drop=FALSE], na.rm=TRUE)
            mono.cc <- is.na(freq) | freq == 0 | freq == 1
            mono <- mono | mono.cc
        }
    }
    mono
}

.CI <- function(Est, SE, CI) {
    LL <- Est + qnorm((1-CI)*0.5)*SE
    UL <- Est + qnorm(1-((1-CI)*0.5))*SE
    c(LL=LL, UL=UL)
}

.waldTest <- function(Est, cov) {
    if (length(Est) == 1) {
        W <- (Est^2)/cov
    } else {
        W <- as.numeric(t(Est) %*% solve(cov) %*% Est)
    }
    pval <- pchisq(W, df=length(Est), lower.tail=FALSE)
    c(Wald.Stat=W, Wald.pval=pval)
}

## give the data frame a crazy name so we don't overwrite something important
.runRegression <- function(model.formula, ..model..data.., model.type, CI, robust, LRtest) {
    tryCatch({
        if (model.type == "linear") {
            mod <- lm(model.formula, data=..model..data..)
        } else if (model.type == "logistic") {
            mod <- glm(model.formula, data=..model..data.., family=binomial())
        } else if (model.type == "poisson") {
            mod <- glm(model.formula, data=..model..data.., family=poisson())
        }
        
        if (!robust) {
            ## model based variance
            Vhat <- vcov(mod)
        } else {
            ## sandwich variance
            Vhat <- vcovHC(mod, type="HC0")
        }
        Est <- unname(coef(mod)["genotype"])
        cov <- Vhat["genotype","genotype"]
        SE <- sqrt(cov)
        ret <- c(Est=Est, SE=SE, .CI(Est, SE, CI), .waldTest(Est, cov))

        if (LRtest) {
            ## promote the data frame to the global environment so lrtest can find it
            ..model..data.. <<- ..model..data..
            mod2 <- lrtest(mod, "genotype")[2,]
            ret <- c(ret, LR.Stat=mod2[["Chisq"]], LR.pval=mod2[["Pr(>Chisq)"]])
            rm(..model..data.., envir=.GlobalEnv)
        }

        ## GxE
        GxE.idx <- grep(":genotype", names(coef(mod)), fixed=TRUE)
        if (length(GxE.idx) > 0) {
            test.GxE <- setNames(.waldTest(coef(mod)[GxE.idx], Vhat[GxE.idx,GxE.idx]),
                                 c("GxE.Stat", "GxE.pval"))
            
            Gj.idx <- grep("genotype", names(coef(mod)), fixed=TRUE)
            test.Joint <- setNames(.waldTest(coef(mod)[Gj.idx], Vhat[Gj.idx,Gj.idx]),
                                   c("Joint.Stat", "Joint.pval"))
            ret <- c(ret, test.GxE, test.Joint)
        }
        
        ret
    }, warning=function(w) NA, error=function(e) NA)
}

.runFirth <- function(model.formula, model.data, CI, PPLtest, geno.index=NULL) {
    tryCatch({
        mod <- logistf(model.formula, data=model.data, alpha=(1-CI),
                       pl=PPLtest, plconf=geno.index, dataout=FALSE)
        ind <- which(mod$terms == "genotype")
        if (PPLtest) stopifnot(ind == geno.index)
        Est <- unname(coef(mod)[ind])
        cov <- vcov(mod)[ind,ind]
        LL <- unname(mod$ci.lower[ind])
        UL <- unname(mod$ci.upper[ind])
        pval <- unname(mod$prob[ind])
        Stat <- qchisq(pval, df=1, lower.tail=FALSE)
        
        ret <- c(Est=Est, SE=sqrt(cov), LL=LL, UL=UL)
        if (PPLtest) {
            c(ret, .waldTest(Est, cov), PPL.Stat=Stat, PPL.pval=pval)
        } else {
            c(ret, Wald.Stat=Stat, Wald.pval=pval)
        }
    }, warning=function(w) NA, error=function(e) NA)
}

    
assocTestReg <- function(genoData,
                         outcome,
                         model.type = c("linear", "logistic", "poisson", "firth"),
                         gene.action = c("additive", "dominant", "recessive"),
                         covar = NULL,
                         ivar = NULL,
                         scan.exclude = NULL,
                         CI = 0.95,
                         robust = FALSE,
                         LRtest = FALSE,
                         PPLtest = TRUE,
                         effectAllele = c("minor", "alleleA"),
                         snpStart = NULL,
                         snpEnd = NULL,
                         block.size = 5000,
                         verbose = TRUE) {

    ## check that arguments are valid
    model.type <- match.arg(model.type)
    gene.action <- match.arg(gene.action)
    effectAllele <- match.arg(effectAllele)
    
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
    dat <- .modelData(genoData, outcome, covar, ivar)
    ## identify samples with any missing data
    keep <- keep & complete.cases(dat)
    dat <- dat[keep,]
    if (!is.null(ivar)) ivar <- paste0(ivar, ":genotype")
    model.formula <- as.formula(paste(outcome, "~", paste(c(covar, ivar, "genotype"), collapse="+")))

    ## for firth test - determine index of genotype in model matrix
    if (model.type == "firth") {
        tmp <- cbind(dat, "genotype"=0)
        geno.index <- which(colnames(model.matrix(model.formula, tmp)) == "genotype")
        rm(tmp)
    }

    ## sample size, assuming no missing genotypes
    n <- nrow(dat)
    if (verbose) message("Running analysis with ", n, " Samples")

    ## number of SNPs in the segment
    nsnp.seg <- snpEnd - snpStart + 1
    ## determine number of SNP blocks
    nblocks <- ceiling(nsnp.seg/block.size)

    ## set up results matrix
    nv <- c("snpID","chr","n","MAF","minor.allele","Est","SE","LL","UL","Wald.Stat","Wald.pval")
    if (LRtest & model.type != "firth") nv <- c(nv, "LR.Stat", "LR.pval")
    if (PPLtest & model.type == "firth") nv <- c(nv, "PPL.Stat", "PPL.pval")
    if (!is.null(ivar) & model.type != "firth") nv <- c(nv, "GxE.Stat", "GxE.pval", "Joint.Stat", "Joint.pval")
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
        geno <- getGenotype(genoData, snp=c(snp.start.pos, nsnp.block), scan=c(1,-1), drop=FALSE)
        ## subset
        geno <- geno[,keep,drop=FALSE]
        
        ## allele frequency
        freq <- .freqFromGeno(genoData, geno, chr[bidx], keep)
        
        ## MAF
        major <- freq > 0.5
        maf <- ifelse(major, 1-freq, freq)
        res[bidx,"MAF"] <- maf
        ## minor allele coding:  A = 1, B = 0
        res[bidx,"minor.allele"] <- ifelse(major, 0, 1)

        ## effect allele
        if (effectAllele == "minor") geno[major,] <- 2 - geno[major,]

        ## transform genotype for gene action
        geno <- .transformGenotype(geno, gene.action)
        
        ## check for monomorphic SNPs
        mono <- .monomorphic(geno, dat[[outcome]], model.type)
        
        ## sample size
        n <- rowSums(!is.na(geno))
        res[bidx, "n"] <- n

        ## loop through SNPs in block
        midx <- (1:nsnp.block)[!mono]
        for (i in midx) {
            mdat <- cbind(dat, genotype=geno[i,])
            mdat <- mdat[complete.cases(mdat),]
            if (model.type %in% c("linear", "logistic")) {
                tmp <- .runRegression(model.formula, mdat, model.type, CI, robust, LRtest)
            } else if (model.type == "firth") {
                tmp <- .runFirth(model.formula, mdat, CI, PPLtest, geno.index)
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
