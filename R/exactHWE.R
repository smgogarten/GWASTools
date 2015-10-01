
## if chr=X, returns keep vector with females only
.checkSexChrHWE <- function(genoData, chr, keep) {
    if (YchromCode(genoData) %in% chr) stop("HWE test not valid for Y chromosome")
    if (MchromCode(genoData) %in% chr) stop("HWE test not valid for mitochondrial SNPs")
    if (XchromCode(genoData) %in% chr) {
        ## check for sex variable
        if (!hasSex(genoData)) {
            stop("Sex values for the samples are required for X chromosome SNPs")
        }
        if (!all(chr == XchromCode(genoData))) {
            stop("X chromosome must be analyzed separately")
        }
        ## only keep females
        keep <- keep & (getSex(genoData) == "F")
    }
    keep
}

  
.countGenotypes <- function(genotypes) {
    nAA <- rowSums(genotypes == 2, na.rm=TRUE)
    nAa <- rowSums(genotypes == 1, na.rm=TRUE)
    naa <- rowSums(genotypes == 0, na.rm=TRUE)
    cbind(nAA, nAa, naa)
}

.permuteGenotypesSlow <- function(genotypes){
  
  counts <- .countGenotypes(genotypes)
  missing <- is.na(genotypes)
  counts <- cbind(counts,
                  nA=2*counts[, "nAA"] + counts[, "nAa"],
                  na=2*counts[, "naa"] + counts[, "nAa"],
                  nmiss=2*rowSums(missing))
  
  
  # slow for loop:
  permuted <- matrix(NA, ncol=ncol(genotypes), nrow=nrow(genotypes))
  for (i in 1:nrow(counts)){
    alleles <- sample(c(rep(1, counts[i, "nA"]), rep(0, counts[i, "na"])))
    genoPermuted <- alleles[c(TRUE, FALSE)] + alleles[c(FALSE, TRUE)]
    permuted[i, !missing[i, ]] <- genoPermuted
  }
  
  permuted

}

.permuteGenotypes <- function(genotypes){
  
  nsamp <- ncol(genotypes)

  alleles <- cbind(genotypes > 0, genotypes > 1)
  shuffled <- t(alleles)
  tmp <- which(!is.na(shuffled), arr.ind=TRUE)
  index <- which(!is.na(shuffled))
  tmp2 <- index[order(tmp[, "col"], runif(nrow(tmp)))]
  # preserve NAs
  ind <- tmp2[match(1:length(shuffled), index)]
  shuffled <- matrix(shuffled[ind], nrow=nrow(alleles), ncol=ncol(alleles), byrow=TRUE) 
  
  permuted <- shuffled[, 1:nsamp] + shuffled[, (nsamp+1):(2*nsamp)]
  permuted
}


exactHWE <- function(genoData,
                     scan.exclude = NULL,
                     geno.counts = TRUE,
                     snpStart = NULL,
                     snpEnd = NULL,
                     block.size = 5000,                      
                     verbose = TRUE,
                     permute=FALSE) {

    ## set snpStart and snpEnd
    if (is.null(snpStart)) snpStart <- 1
    if (is.null(snpEnd)) snpEnd <- nsnp(genoData)

    ## set which samples to keep
    keep <- .keepSamples(genoData, scan.exclude)

    ## get chromosome information
    chr <- getChromosome(genoData, index=snpStart:snpEnd)

    ## sex chromosome checks
    keep <- .checkSexChrHWE(genoData, chr, keep)

    ## number of SNPs in the segment
    nsnp.seg <- snpEnd - snpStart + 1
    nblocks <- ceiling(nsnp.seg/block.size)

    ## set up results matrix
    nv <- c("snpID", "chr")
    if (geno.counts) nv <- c(nv, "nAA", "nAB", "nBB")
    nv <- c(nv, "MAF", "minor.allele", "f", "pval")
    res <- matrix(NA, nrow=nsnp.seg, ncol=length(nv), dimnames=list(NULL, nv))
    
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

        if (permute) geno <- .permuteGenotypes(geno)
        
        ## count genotypes
        tmpGenotypeCounts <- .countGenotypes(geno)
        if (geno.counts) {
            res[bidx,c("nAA", "nAB", "nBB")] <- tmpGenotypeCounts
        }

        ## allele frequency
        freq <- 0.5*rowMeans(geno, na.rm=TRUE)
        major <- freq > 0.5 & !is.na(freq)
        maf <- ifelse(major, 1-freq, freq)
        res[bidx,"MAF"] <- maf
        ## minor allele coding:  A = 1, B = 0
        res[bidx,"minor.allele"] <- ifelse(major, 0, 1)
        
        ## calculate inbreeding coefficient
        obs.het <- tmpGenotypeCounts[,2]
        geno.tot <- rowSums(tmpGenotypeCounts)
        exp.het <- 2*maf*(1-maf)*geno.tot
        res[bidx,"f"] <- 1-(obs.het/exp.het)
      
        ## calculate HW p-val
        res[bidx,"pval"] <- HWExact(as.data.frame(tmpGenotypeCounts))
        
        rate <- format(Sys.time() - startTime, digits=4)        
        if (verbose) message(paste("Block", b, "of", nblocks, "Completed -", rate))
    }

    ## results data frame
    res <- as.data.frame(res)
    res$snpID <- getSnpID(genoData, index=snpStart:snpEnd)

    ## convert minor.allele coding back to A/B
    res[,"minor.allele"][res[,"minor.allele"] == 1] <- "A"
    res[,"minor.allele"][res[,"minor.allele"] == 0] <- "B"

    ## set pvalue to NA for monomorphic SNPs
    res$pval[res$MAF == 0] <- NA
  
    return(res)
}
