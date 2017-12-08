
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
            write.table(out, con, sep="\t", col.names=FALSE, row.names=FALSE,
                        quote=FALSE)
        }
    }

    ## close output file
    close(con)
}



vcfCheck <- function(genoData, vcf.file, sample.col="scanID", id.col="snpID",
                     scan.exclude=NULL, snp.exclude=NULL,
                     block.size=1000, verbose=TRUE) {

    # check for samp identifier
    stopifnot(hasScanVariable(genoData, sample.col))
    # make 2 col df with samp identifier and scanID (may be =)
    samp.id.dat <- getScanVariable(genoData, c(sample.col, "scanID"))

    # check for snp identifier 
    stopifnot(hasSnpVariable(genoData, id.col))
    # make 2 col df with snp identifier and snpID (may be =)        
    snp.id.dat <- getSnpVariable(genoData, c(id.col, "snpID"))
    
    alleleA <- getAlleleA(genoData)

    ## subset to expected snpIDs and scanIDs based on exclude args
    ## note we're allowing for VCF to have more data than genoData
    if (!is.null(scan.exclude)) {
      if(verbose){message("excluding ", prettyNum(length(scan.exclude), big.mark=","),
                          " samples from check")}
      scan.include <- setdiff(getScanID(genoData), scan.exclude)
    } else {
      scan.include <- getScanID(genoData)
    }
    
    if (!is.null(snp.exclude)) {
      if(verbose){message("excluding ", prettyNum(length(snp.exclude), big.mark=","),
                          " SNPs from check")}
      snp.include <- setdiff(getSnpID(genoData), snp.exclude)
    } else {
      snp.include <- getSnpID(genoData)
    }

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
    samples <- header[10:ncol]

    # check if there are samples in VCF and not genoData
    # give message for VCF overall, vs. for each geno block below
    vcf.samp.only <- setdiff(samples, samp.id.dat[,sample.col])
    if(length(vcf.samp.only) > 0){
      message("Note VCF has ", prettyNum(length(vcf.samp.only), big.mark=","),
              " sample(s) not present in genoData;\n these will be excluded from the check")
    }

    # reads VCF file block.size lines at a time
    check.total <- 0
    while (length(x <- scan(vcf, what=character(), nlines=block.size, quiet=TRUE)) > 0) {
        geno.vcf <- matrix(x, ncol=ncol, byrow=TRUE)
        id <- geno.vcf[,3]
        ref.vcf <- geno.vcf[,4]
        geno.vcf <- geno.vcf[,10:ncol]

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
        dimnames(geno.vcf) <- list(id, samples)

        #### allow VCF with more data than present in genoData
        ## subset on SNPs present in genoData
        vcf.snp.only <- setdiff(id, snp.id.dat[,id.col])

        ## create snp id vector for subsetting VCF
        if(length(vcf.snp.only) > 0){
          message("Note VCF block has ", prettyNum(length(vcf.snp.only), big.mark=","),
              " SNP(s) not present in genoData;\n these will be excluded from the check")
          ref.vcf <- ref.vcf[id %in% snp.id.dat[,id.col]]
          id <- id[id %in% snp.id.dat[,id.col]]          
        }

        ## create samp id vector for subsetting VCF
        ## extra samps were already reported before while loop
        ## do this silently
        if(length(vcf.samp.only) > 0) {
          samples <- samples[samples %in% samp.id.dat[,sample.col]]
        }

        # if 'id' or 'samples' lists were not reduced, then no subsetting occurs
        geno.vcf <- geno.vcf[id, samples]
        
        ## add check for expected snps and samples

        ## check is samples are in same order - give warning in diff order

        ## match order before doing comparison
        
        ## subset genoData to requested snps, samples
        ## as specified by exclude arguments

        count <- nrow(geno.vcf)
        start.orig <- which(snp.id == id[1])
        end.orig <- which(snp.id == id[count])
        count.orig <- end.orig - start.orig + 1
        geno.orig <- getGenotype(genoData, scan=c(1,-1), snp=c(start.orig,count.orig))
        dimnames(geno.orig) <- list(snp.id[start.orig:end.orig], samp.id)
        geno.orig <- geno.orig[rownames(geno.vcf), colnames(geno.vcf)]
        ref.orig <- alleleA[match(rownames(geno.vcf), snp.id)]
        allele.switch <- ref.orig != ref.vcf
        geno.orig[allele.switch,] <- 2 - geno.orig[allele.switch,]

        # matrix of all genos - made to be in same order snp and sample-wise
        stopifnot(allequal(geno.vcf, geno.orig))

        check.total <- check.total + count
        message("Checked ", check.total, " SNPs")
    }

    close(vcf)
}
