## 1. To use an GxE interaction term, you MUST specify the E variable with covars AND ivar.
##    This is because the notation used is G:E which does not add E to the model, whereas G*E does add G to the model
##    G is assumed by default with any interaction term
##
## 2. The code below is set up so that "mod" is the full model, including the genotype, covariates and the
##    interaction term if one is specified and (a) mod2 is a reduced model with no genotype term in order
##    to calculate the LR p-value for the genotype and (b) mod3 is a reduced model with no interaction term
##    in order to calculate the LR p-value for the interaction term
##
## 3. If an interaction term is specified, the Beta and standard error are retrieved from mod which is the full model

.coxLR <- function(surv, mod, model.string, ..model.data..) {
    # if testing genotype only return overall likelihood ratio computed by coxph
    if(length(coef(mod))==1) {
        sumcph  <- summary(mod)
        LR.Stat <- sumcph$logtest[1]
        LR.pval <- 1 - pchisq(LR.Stat, 1)
    } else {
        # if covariates, create reduced model without genotype and compare with 
        terms  <- unlist(strsplit(model.string, " + ", fixed=TRUE))
        model2 <- as.formula(paste(terms[!grepl("genotype", terms, fixed=TRUE)], collapse=" + "))
        # be sure to remove scans that have missing genotypes from the reduced model
        # so can properly compare with same samples sizes in anova() call below
        genotype <- NULL # hack to use subset arg below
        mod2   <- coxph(model2, data=..model.data.., subset=(!is.na(genotype)))
        resLR  <- anova(mod,mod2)               
        LR.Stat <- resLR[["Chisq"]][2]
        # name of column changed in R 4.2.1
        p.col.names <- c("P(>|Chi|)", "Pr(>|Chi|)")
        p.col <- intersect(names(resLR), p.col.names)
        LR.pval <- resLR[[p.col]][2]
        ## explicit computation code if even needed
        ## LR.Stat <- -2*(mod2$loglik[2] - mod$loglik[2])
        ## LR.pval <- 1-pchisq(LR.Stat,1)     
    }
    c(LR.Stat=LR.Stat, LR.pval=LR.pval)
}

## GxE
# use index of interaction term from "mod", the full model, to (a) drop the interaction term in "mod3" in order to
# calculate the LR p-value of the GxE term and (b) retrieve Est and SE for the interaction term from "mod"
.coxGxE <- function(surv, mod, model.string, ..model.data..) {
    GxE.idx <- grep(":genotype", names(coef(mod)), fixed=TRUE)
    if (length(GxE.idx) == 0) return(NULL)
    terms   <- unlist(strsplit(model.string, " + ", fixed=TRUE))
    model3  <- as.formula(paste(terms[!grepl(":genotype", terms, fixed=TRUE)], collapse=" + "))
    mod3    <- coxph(model3, data=..model.data..)
    lr      <- -2*(mod3$loglik[2] - mod$loglik[2])
    pval    <- 1-pchisq(lr, 1)
    if (length(GxE.idx) == 1) {
        GxE.Est <- unname(coef(mod)[GxE.idx])
        GxE.SE  <- sqrt(vcov(mod)[GxE.idx,GxE.idx])
        ret <- c(GxE.Est=GxE.Est, GxE.SE=GxE.SE, GxE.Stat=lr, GxE.pval=pval)
    } else {
        ret <- c(GxE.Est=NA, GxE.SE=NA, GxE.Stat=lr, GxE.pval=pval)
    }
    ret
}

.runCPH <- function(model.string, ..model.data.., event, time.to.event, LRtest=FALSE) {
    tryCatch(
    {
        surv <- Surv(time=..model.data..[[time.to.event]], event=..model.data..[[event]])
        model.formula <- as.formula(model.string)
        mod  <- coxph(model.formula, data=..model.data..)
        Est  <- unname(coef(mod)["genotype"])
        cov  <- vcov(mod)["genotype","genotype"]
        SE   <- sqrt(vcov(mod)["genotype","genotype"])
        resWald <- .waldTest(Est, cov)
        ret <- c(Est=Est, SE=SE, resWald[c("Wald.Stat", "Wald.pval")])

        if (LRtest) {
            lrtest <- .coxLR(surv, mod, model.string, ..model.data..)
            ret <- c(ret, lrtest)
        }
        
        gxe <- .coxGxE(surv, mod, model.string, ..model.data..)
        ret <- c(ret, gxe)
        ret
    }, # matches tryCatch {
    warning=function(w) {
        # should have been sent here from "mod  <- coxph(model.formula, data=..model.data..)" above assume
        # we would not have the case where "mod  <- coxph(" passes, but "mod2  <- coxph(" from LR test fails
        # need to run the test that had the exception again, but this time outside of tryCatch since none of
        # mod, Est, or SE has state in this block
        mod  <- coxph(model.formula, data=..model.data..)
        if (!LRtest) {
            # we assume only LR test will be valid, so if not set return NA 
            ret <- NA
        } else {
            print("hi") 
            Est  <- unname(coef(mod)["genotype"])
            SE   <- sqrt(vcov(mod)["genotype","genotype"])
            lrtest <- .coxLR(surv, mod, model.string, ..model.data..)
            print(lrtest)
            ret <- c(Est=Est, SE=SE, Wald.Stat=NA, Wald.pval=NA, lrtest)

            gxe <- .coxGxE(surv, mod, model.string, ..model.data..)
            ret <- c(ret, gxe)
            ret
        }
        
    }, # matches warning=function(w) {
    error=function(e) NA
    ) # matches tryCatch (
}

assocCoxPH <- function(genoData,
                       event,
                       time.to.event,
                       gene.action = c("additive", "dominant", "recessive"),
                       covar = NULL,
                       ivar = NULL,
                       strata = NULL,
                       cluster = NULL,
                       scan.exclude = NULL,
                       LRtest = FALSE,
                       effectAllele = c("minor", "alleleA"),
                       snpStart = NULL,
                       snpEnd = NULL,
                       block.size = 5000,
                       verbose = TRUE) {

    ## check that arguments are valid
    gene.action <- match.arg(gene.action)
    effectAllele <- match.arg(effectAllele)
    
    ## set snpStart and snpEnd
    if (is.null(snpStart)) snpStart <- 1
    if (is.null(snpEnd)) snpEnd <- nsnp(genoData)

    ## set which samples to keep
    keep <- .keepSamples(genoData, scan.exclude)

    ## get chromosome information
    chr <- getChromosome(genoData, index=snpStart:snpEnd)

    ## sex chromosome checks
    keep <- .checkSexChr(genoData, chr, keep)

    ## read in outcome and covariate data
    if (verbose) message("Reading in Phenotype and Covariate Data...")
    dat <- .modelData(genoData, chr, event, c(time.to.event, covar, strata, cluster), ivar)
    ## identify samples with any missing data
    keep <- keep & complete.cases(dat)
    dat <- dat[keep,,drop=FALSE]
    if (!is.null(ivar)) ivar <- paste0(ivar, ":genotype")
    if (!is.null(strata)) strata <- paste0("strata(", paste(strata, collapse=","), ")")
    if (!is.null(cluster)) cluster <- paste0("cluster(", paste(cluster, collapse=","), ")")
    model.string <- paste("surv ~", paste(c(covar, ivar, strata, cluster, "genotype"), collapse=" + "))

    ## sample size, assuming no missing genotypes
    n <- nrow(dat)
    if (verbose) message("Running analysis with ", n, " Samples")

    ## number of SNPs in the segment
    nsnp.seg <- snpEnd - snpStart + 1
    nblocks <- ceiling(nsnp.seg/block.size)

    ## set up results matrix
    nv <- c("snpID", "chr", "n", "n.events", "effect.allele", "EAF", "MAF", "maf.filter",
            "Est", "SE", "Wald.Stat", "Wald.pval")
    if (LRtest) nv <- c(nv, "LR.Stat", "LR.pval")
    if (!is.null(ivar)) nv <- c(nv, "GxE.Est", "GxE.SE", "GxE.Stat", "GxE.pval")
    res <- matrix(NA, nrow=nsnp.seg, ncol=length(nv), dimnames=list(NULL, nv))
    reg.cols <- which(colnames(res) == "Est"):ncol(res)

    ## chromosome
    res[,"chr"] <- chr

    if (verbose) message("Beginning Calculations...")
    for (b in 1:nblocks) {
        
        ## keep track of time for rate reporting
        startTime <- Sys.time()
        
        snp.start.pos <- snpStart + (b-1)*block.size
        nsnp.block <- ifelse(snp.start.pos + block.size > snpEnd,
                             snpEnd - snp.start.pos + 1, block.size)
        bidx <- ((b-1)*block.size + 1):((b-1)*block.size + nsnp.block)
        
        ## get genotypes for the block
        geno <- getGenotype(genoData, snp=c(snp.start.pos, nsnp.block), scan=c(1,-1), drop=FALSE)
        geno <- geno[,keep,drop=FALSE]
        
        ## allele frequency
        freq <- .freqFromGeno(genoData, geno, chr[bidx], keep)
        major <- freq > 0.5 & !is.na(freq)
        maf <- ifelse(major, 1-freq, freq)
        res[bidx,"MAF"] <- maf

        ## effect allele
        if (effectAllele == "minor") {
            geno[major,] <- 2 - geno[major,]
            ## minor allele coding:  A = 1, B = 0
            res[bidx,"effect.allele"] <- ifelse(major, 0, 1)
            res[bidx,"EAF"] <- maf
        } else {
            res[bidx,"effect.allele"] <- 1
            res[bidx,"EAF"] <- freq
        }

        ## transform genotype for gene action
        geno <- .transformGenotype(geno, gene.action)
        
        ## check for monomorphic SNPs
        mono <- .monomorphic(geno, dat[[event]], "survival")

        ## sample size
        res[bidx, "n"] <- rowSums(!is.na(geno))
        ne <- rowSums(!is.na(geno[,as.logical(dat[[event]])]))
        res[bidx, "n.events"] <- ne
        
        ## filter
        ## calculate MAF with male dosage=2, even for X chrom
        maf2 <- 0.5*rowMeans(geno, na.rm=TRUE)
        maf2 <- pmin(maf2, 1-maf2)
        res[bidx, "maf.filter"] <- 2*maf2*(1-maf2)*ne > 75
        
        ## loop through SNPs in block
        midx <- (1:nsnp.block)[!mono]
        for (i in midx) {
            mdat <- cbind(dat, genotype=geno[i,])
            mdat <- mdat[complete.cases(mdat),]            
            tmp <- .runCPH(model.string, mdat, event, time.to.event, LRtest)
            res[bidx[i], reg.cols] <- tmp
        }
        
        rate <- format(Sys.time() - startTime, digits=4)        
        if (verbose) message(paste("Block", b, "of", nblocks, "Completed -", rate))
    }

    ## results data frame
    res <- as.data.frame(res)
    res$snpID <- getSnpID(genoData, index=snpStart:snpEnd)
    res$maf.filter <- as.logical(res$maf.filter)

    ## convert effect.allele coding back to A/B
    res$effect.allele[res$effect.allele == 1] <- "A"
    res$effect.allele[res$effect.allele == 0] <- "B"

    return(res)
}
