
## should replace this with a method!
.getMetadata <- function(object, varname) {
    varMetadata(object@snpAnnot)[varname, "labelDescription"]
}

.setFilter <- function(genoData, filter.cols) {
    filt.df <- data.frame(getSnpVariable(genoData, filter.cols))
    names(filt.df) <- filter.cols
    filter <- rep("", nsnp(genoData))
    meta <- character()
    for (c in filter.cols) {
        filter[filt.df[[c]]] <- paste(filter[filt.df[[c]]], c, sep=";")
        meta[c] <- paste0('##FILTER=<ID=', c,
                          ',Description="', .getMetadata(genoData, c), '">')
    }
    filter <- sub("^;", "", filter)
    filter[filter == ""] <- "PASS"
    list(filter=filter, meta=meta)
}

.setInfo <- function(genoData, info.cols) {
    info.df <- data.frame(getSnpVariable(genoData, info.cols))
    names(info.df) <- info.cols
    info <- rep("", nsnp(genoData))
    meta <- character()
    for (c in info.cols) {
        if (is.logical(info.df[[c]])) {
            info[info.df[[c]]] <- paste(info[info.df[[c]]], c, sep=";")
        } else {
            info <- paste(info, paste0(c, "=", info.df[[c]]), sep=";")
        }
        type.map <- c(integer="Integer", numeric="Float", character="String",
                      factor="String", logical="Flag")
        meta[c] <- paste0('##INFO=<ID=', c,
                          ',Number=', ifelse(is.logical(info.df[[c]]), 0, 1),
                          ',Type=', type.map[class(info.df[[c]])],
                          ',Description="', .getMetadata(genoData, c), '">')
    }
    info <- sub("^;", "", info)
    info[info == ""] <- "."
    list(info=info, meta=meta)
}

vcfWrite <- function(genoData, vcf.file="out.vcf", sample.col="scanID",
                     id.col="snpID", qual.col=NULL, filter.cols=NULL,
                     info.cols=NULL, scan.exclude=NULL, snp.exclude=NULL,
                     scan.order=NULL,
                     ref.allele=NULL, block.size=1000, verbose=TRUE) {
    ## fixed fields
    chrom <- getChromosome(genoData, char=TRUE)
    pos <- getPosition(genoData)
    id <- getSnpVariable(genoData, id.col)
    ## check for missing values in id
    id <- as.character(id)
    id[is.na(id) | id == ""] <- "."
    if (!is.null(ref.allele)) {
        stopifnot(length(ref.allele) == nsnp(genoData))
        stopifnot(all(ref.allele %in% c("A", "B")))
        a <- getAlleleA(genoData)
        b <- getAlleleB(genoData)
        ref <- ifelse(ref.allele == "A", a, b)
        alt <- ifelse(ref.allele == "A", b, a)
    } else {
        ref <- getAlleleA(genoData)
        alt <- getAlleleB(genoData)
    }
    if (is.null(qual.col)) {
        qual <- rep(".", nsnp(genoData))
    } else {
        qual <- getSnpVariable(genoData, qual.col)
    }
    if (is.null(filter.cols)) {
        filter <- rep(".", nsnp(genoData))
        filt.meta <- character()
    } else {
        filt.list <- .setFilter(genoData, filter.cols)
        filter <- filt.list[["filter"]]
        filt.meta <- filt.list[["meta"]]
        rm(filt.list)
    }
    if (is.null(info.cols)) {
        info <- rep(".", nsnp(genoData))
        info.meta <- character()
    } else {
        info.list <- .setInfo(genoData, info.cols)
        info <- info.list[["info"]]
        info.meta <- info.list[["meta"]]
        rm(info.list)
    }
    format <- rep("GT", nsnp(genoData))
    fixed <- cbind(chrom, pos=as.character(pos), id=id,
                   ref, alt, qual=as.character(qual), filter,
                   info, format)
    ## subset with snp.exclude
    if (!is.null(snp.exclude)) {
        snp.index <- !(getSnpID(genoData) %in% snp.exclude)
    } else {
        snp.index <- rep(TRUE, nsnp(genoData))
    }

    ## open output file
    con <- file(vcf.file, "w")

    ## write metadata
    meta <- c('##fileformat=VCFv4.1',
              paste0('##fileDate=', Sys.Date()),
              '##source=GWASTools::vcfWrite()',
              filt.meta,
              info.meta,
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    writeLines(meta, con)

    ## get samples
    samples <- getScanVariable(genoData, sample.col)
    scanID <- getScanID(genoData)
    if (is.null(scan.order)) {
        scan.order <- scanID
    } else {
        stopifnot(all(scan.order %in% scanID))
    }
    ## subset with scan.exclude
    if (!is.null(scan.exclude)) {
        scan.order <- setdiff(scan.order, scan.exclude)
    }
    ## put samples in new order if requested
    scan.index <- match(scan.order, scanID)

    ## write header
    hdr <- paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                   samples[scan.index]), collapse="\t")
    writeLines(hdr, con)
    
    ## fwrite cannot write to an open connection
    close(con)

    ## write genotypes in blocks
    nblocks <- ceiling(nsnp(genoData) / block.size)
    for (i in 1:nblocks) {
        start <- (i-1)*block.size + 1
        end <- min(i*(block.size), nsnp(genoData))
        count <- end - start + 1
        n <- sum(snp.index[start:end])
        if (verbose) message("Block ", i, " of ", nblocks, "... ", n, " SNPs")

        if (n > 0) {
            geno <- getGenotype(genoData, snp=c(start,count), scan=c(1,-1))
            ## switch allele coding if ref.allele is B
            geno[ref.allele[start:end] == "B",] <- 2 - geno[ref.allele[start:end] == "B",]
            geno <- geno[snp.index[start:end], scan.index, drop=FALSE]
            geno[is.na(geno)] <- "./."
            geno[geno == 2] <- "0/0"
            geno[geno == 1] <- "0/1"
            geno[geno == 0] <- "1/1"

            out <- cbind(fixed[(start:end)[snp.index[start:end]],,drop=FALSE], geno)
            fwrite(out, vcf.file, append=TRUE, quote=FALSE, sep="\t",
                   row.names=FALSE, col.names=FALSE)
        }
    }
}


vcfCheck <- function(genoData, vcf.file, sample.col="scanID", id.col="snpID",
                     scan.exclude=NULL, snp.exclude=NULL,
                     block.size=1000, verbose=TRUE) {

    ## check for samp identifier
    stopifnot(hasScanVariable(genoData, sample.col))
    ## make 2 col df with samp identifier and scanID (may be equivalent)
    samp.id.dat <- getScanVariable(genoData, c(sample.col, "scanID"))

    ## check for snp identifier 
    stopifnot(hasSnpVariable(genoData, id.col))
    ## make 2 col df with snp identifier and snpID (may be equivalent)        
    snp.id.dat <- getSnpVariable(genoData, c(id.col, "snpID"))
    
    alleleA <- getAlleleA(genoData)

    ## get vectors of expected snpIDs and scanIDs based on exclude args
    ## note we're allowing for VCF to have more data than genoData
    if (!is.null(scan.exclude)) {
      message("Excluding ", prettyNum(length(scan.exclude), big.mark=","),
              " genoData samples from check")
      scan.include <- setdiff(getScanID(genoData), scan.exclude)
    } else {
      scan.include <- getScanID(genoData)
    }
    ## add logical flag to samp.id.dat to indicate inclusion in check
    samp.id.dat$include <- samp.id.dat$scanID %in% scan.include
    
    if (!is.null(snp.exclude)) {
      message("Excluding ", prettyNum(length(snp.exclude), big.mark=","),
               " genoData SNPs from check")
      snp.include <- setdiff(getSnpID(genoData), snp.exclude)
    } else {
      snp.include <- getSnpID(genoData)
    }
    ## add logical flag to snp.id.dat to indicate inclusion in check
    snp.id.dat$include <- snp.id.dat$snpID %in% snp.include    

    ## open VCF
    vcf <- file(vcf.file, "r")

    ## read until we get the header
    while (length(s <- readLines(vcf, n=1)) > 0) {
        if (substr(s, 1, 6) == "#CHROM") {
            header <- scan(text=s, what=character(), quiet=TRUE)
            break
        }
    }
    ncol <- length(header)
    samp.allvcf <- header[10:ncol]

    ## relay if there are samples in VCF and not genoData
    vcf.samp.only <- setdiff(samp.allvcf, samp.id.dat[,sample.col])
    if(length(vcf.samp.only) > 0){
      message("Note VCF has ", prettyNum(length(vcf.samp.only), big.mark=","),
              " sample(s) not present in genoData;\nthese will be excluded from the check")
    }

    ## check if VCF includes more genoData samps than are retained
    ## after applying exclude lists
    ## ...given common scenario where vcfCheck is run following vcfWrite,
    ## and the two functions are given the same exclude lists
    samp.id.dat$vcf <- samp.id.dat[, sample.col] %in% samp.allvcf
    if(sum(with(samp.id.dat, vcf & !include)) > 0){
      message("Note after applying exclude argument, VCF contains additional genoData samples")
    }

    ## write vector of VCF sample IDs to include in check
    ## in VCF file order
    samp.vcf.incl <- samp.allvcf[samp.allvcf %in%
                                 with(samp.id.dat, samp.id.dat[vcf & include, sample.col])]
    
    ## check if samples are in same order - give warning in diff order
    if(!allequal(as.character(samp.vcf.incl),
                 as.character(samp.id.dat[samp.id.dat$include, sample.col]))) {
      message("Note, sample order in VCF differs from genoData")
    }

    ## reads VCF file block.size lines at a time
    check.total <- 0
    while (length(x <- scan(vcf, what=character(), nlines=block.size, quiet=TRUE)) > 0) {
        geno.vcf <- matrix(x, ncol=ncol, byrow=TRUE)
        id <- geno.vcf[,3]
        ref.vcf <- geno.vcf[,4]
        geno.vcf <- geno.vcf[, 10:ncol, drop=FALSE]

        ## create df of VCF snp ID and ref allele to facilitate
        ## subsequent subsetting of ref allele to match checked SNPs
        ref.dat <- data.frame(var=id, ref=ref.vcf)

        ## take the first 3 characters - GT field for diploid genotypes
        geno.vcf <- substr(geno.vcf, 1, 3)
        ## if phased, convert to unphased separator
        geno.vcf <- sub("|", "/", geno.vcf, fixed=TRUE)
        ## convert to count of REF
        geno.vcf[geno.vcf == "0/0"] <- 2
        geno.vcf[geno.vcf == "0/1"] <- 1
        geno.vcf[geno.vcf == "1/0"] <- 1
        geno.vcf[geno.vcf == "1/1"] <- 0
        geno.vcf[geno.vcf == "./."] <- NA
        mode(geno.vcf) <- "integer"
        dimnames(geno.vcf) <- list(id, samp.allvcf)

        #### allow VCF with more data than present in genoData (prior to exclusions)
        ## give message if running verbose=T
        vcf.snp.only <- setdiff(id, snp.id.dat[,id.col])
        if(length(vcf.snp.only) > 0 & verbose) {
            message("Note VCF block has ", prettyNum(length(vcf.snp.only), big.mark=","),
              " SNP(s) not present in genoData;\nthese will be excluded from the check")
        }

        ## get vector of snps to check
        snp.block.incl <- id[id %in% snp.id.dat[snp.id.dat$include, id.col]]

        ## if verbose, give message if VCF includes more genoData snps 
        ## than are retained after applying exclude lists
        ## (note sample-level check happened before iterating over blocks)
        if(length(setdiff(id, snp.block.incl)) > 0 & verbose) {
           message("Note after applying exclude argument, VCF block contains additional genoData SNPs")
        }

        ## subset block to snps and samples to be checked, i.e.
        ## those overlapping with genoData and not in exclude lists
        geno.vcf <- geno.vcf[snp.block.incl, samp.vcf.incl, drop=FALSE]

        ## only continue if there are snps left to check after exclusions
        if(nrow(geno.vcf) > 0) {

            ## get list of snpIDs to check in this block
            ## could this be done by matching snp.id.dat to geno.vcf rownames?
            ## yes, but coerce to character so it's not used as index
            snp.block.sel <- snp.id.dat[, id.col] %in% snp.block.incl
            snpID.block.incl <- snp.id.dat[snp.block.sel, "snpID"]

            ## get VCF ref allele for snps in this block
            ref.block <- ref.dat$ref[match(snp.block.incl, ref.dat$var)]

            ## using getGenotypeSelection, with snpID and scanID values as args        
            geno.orig <- getGenotypeSelection(genoData,
                                              scanID=samp.id.dat$scanID[samp.id.dat$include],
                                              snpID=snpID.block.incl)

            ## compare alleleA and identify switches (i.e., diff allele being counted)
            ref.orig <- alleleA[match(rownames(geno.vcf), snp.id.dat[, id.col])]
            allele.switch <- ref.orig != ref.block
            if(sum(allele.switch) > 0) {
              ## convert to count same allele, where alleleA/REF differs
              geno.orig[allele.switch,] <- 2 - geno.orig[allele.switch,]
             }

            ## now we have matrices of all genos,
            ## made to be in same order snp and sample-wise
            ## check equality - exit w/informative error if not equal
            ## note 'allequal' checks if NA placement matches
            if(!allequal(geno.vcf, geno.orig)) {

              ## to check if na's are equal, replace NA with -1
              geno.orig[is.na(geno.orig)] <- -1
              geno.vcf[is.na(geno.vcf)] <- -1

              ## this returns logical matrix, with same dimnames as (geno.vcf)
              eq.matx <- geno.vcf == geno.orig

              ## find row index where discrepancies start
              err.row <- which(!eq.matx, arr.ind=TRUE)[1]

              ## report variant name
              err.var <- rownames(geno.vcf)[err.row]
              
              # give
              stop(paste("genoData and VCF are not equal, starting at VCF variant ID",
                         err.var))
            }

            count <- nrow(geno.vcf)
            check.total <- check.total + count
            if(verbose) {message("Checked ", check.total, " SNPs")}
      }
    }

    close(vcf)
}
