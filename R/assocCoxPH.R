.runCPH <- function(model.string, model.data, event, time.to.event) {
    tryCatch({
        surv <- Surv(time=model.data[[time.to.event]], event=model.data[[event]])
        model.formula <- as.formula(model.string)
        mod <- coxph(model.formula, data=model.data)
        Est <- unname(coef(mod)["genotype"])
        SE <- sqrt(vcov(mod)["genotype","genotype"])
        z <- Est/SE
        pval <- 1 - pchisq(z^2, 1)
        ret <- c(Est=Est, SE=SE, z.Stat=z, z.pval=pval)
        
        ## GxE
        GxE.idx <- grep(":genotype", names(coef(mod)), fixed=TRUE)
        if (length(GxE.idx) > 0) {
            terms <- unlist(strsplit(model.string, " + ", fixed=TRUE))
            model2 <- as.formula(paste(terms[!grepl(":genotype", terms, fixed=TRUE)], collapse=" + "))
            mod2 <- coxph(model2, data=model.data)
            lr <- -2*(mod2$loglik[2] - mod$loglik[2])
            pval <- 1-pchisq(lr, 1)
            ret <- c(ret, GxE.Stat=lr, GxE.pval=pval)
        }
        ret
    }, warning=function(w) NA, error=function(e) NA)
}


assocCoxPH <- function(genoData,
                       event,
                       time.to.event,
                       gene.action = c("additive", "dominant", "recessive"),
                       covar = NULL,
                       ivar = NULL,
                       strata = NULL,
                       scan.exclude = NULL,
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
    dat <- .modelData(genoData, chr, event, c(time.to.event, covar, strata), ivar)
    ## identify samples with any missing data
    keep <- keep & complete.cases(dat)
    dat <- dat[keep,]
    if (!is.null(ivar)) ivar <- paste0(ivar, ":genotype")
    if (!is.null(strata)) strata <- paste0("strata(", paste(strata, collapse=","), ")")
    model.string <- paste("surv ~", paste(c(covar, ivar, strata, "genotype"), collapse=" + "))


    ## sample size, assuming no missing genotypes
    n <- nrow(dat)
    if (verbose) message("Running analysis with ", n, " Samples")

    ## number of SNPs in the segment
    nsnp.seg <- snpEnd - snpStart + 1
    nblocks <- ceiling(nsnp.seg/block.size)

    ## set up results matrix
    nv <- c("snpID", "chr", "n.events", "effect.allele", "EAF", "MAF", "filter",
            "Est", "SE", "z.Stat", "z.pval")
    if (!is.null(ivar)) nv <- c(nv, "GxE.Stat", "GxE.pval")
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
        mono <- .monomorphic(geno, dat[[event]], "logistic")

        ## sample size
        ne <- rowSums(!is.na(geno[,as.logical(dat[[event]])]))
        res[bidx, "n.events"] <- ne
        
        ## filter
        ## calculate MAF with male dosage=2, even for X chrom
        maf2 <- 0.5*rowMeans(geno, na.rm=TRUE)
        maf2 <- pmin(maf2, 1-maf2)
        res[bidx, "filter"] <- 2*maf2*(1-maf2)*ne > 75
        
        ## loop through SNPs in block
        midx <- (1:nsnp.block)[!mono]
        for (i in midx) {
            mdat <- cbind(dat, genotype=geno[i,])
            mdat <- mdat[complete.cases(mdat),]            
            tmp <- .runCPH(model.string, mdat, event, time.to.event)
            res[bidx[i], reg.cols] <- tmp
        }
        
        rate <- format(Sys.time() - startTime, digits=4)        
        if (verbose) message(paste("Block", b, "of", nblocks, "Completed -", rate))
    }

    ## results data frame
    res <- as.data.frame(res)
    res$snpID <- getSnpID(genoData, index=snpStart:snpEnd)
    res$filter <- as.logical(res$filter)

    ## convert effect.allele coding back to A/B
    res$effect.allele[res$effect.allele == 1] <- "A"
    res$effect.allele[res$effect.allele == 0] <- "B"

    return(res)
}
