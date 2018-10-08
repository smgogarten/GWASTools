.makeFilter <- function(x, filter.cols) {
    filt.df <- data.frame(getSnpVariable(x, filter.cols))
    names(filt.df) <- filter.cols
    filter <- rep("", nsnp(x))
    meta <- character()
    for (c in filter.cols) {
        filter[filt.df[[c]]] <- paste(filter[filt.df[[c]]], c, sep=";")
    }
    meta <- do.call(rbind,
                    lapply(filter.cols, function(c) {
                        S4Vectors::DataFrame(Description=.getMetadata(x, c), row.names=c)
                    }))
    filter <- sub("^;", "", filter)
    filter[filter == ""] <- "PASS"
    attr(filter, "header") <- meta
    filter
}

.refalt <- function(x, ref.allele) {
    if (!is.null(ref.allele)) {
        stopifnot(length(ref.allele) == nsnp(x))
        stopifnot(all(ref.allele %in% c("A", "B")))
        a <- getAlleleA(x)
        b <- getAlleleB(x)
        ref <- ifelse(ref.allele == "A", a, b)
        alt <- ifelse(ref.allele == "A", b, a)
    } else {
        ref <- getAlleleA(x)
        alt <- getAlleleB(x)
    }
    list(ref=ref, alt=alt)
}

.makeRowRanges <- function(x, id.col, ref.allele) {
    ## fixed fields
    chrom <- getChromosome(x, char=TRUE)
    pos <- getPosition(x)
    id <- getSnpVariable(x, id.col)
    ## check for missing values in id
    id <- as.character(id)
    missing <- is.na(id) | id %in% c("", ".")
    alleles <- .refalt(x, ref.allele)
    id[missing] <- sprintf("%s:%d_%s/%s", chrom[missing], pos[missing],
                           alleles$ref[missing], alleles$alt[missing])
    GenomicRanges::GRanges(chrom, IRanges::IRanges(pos, width=1L, names=id))
}

.makeFixed <- function(x, ref.allele, qual.col, filter.cols) {
    alleles <- .refalt(x, ref.allele)
    ref <- Biostrings::DNAStringSet(alleles$ref)
    alt <- Biostrings::DNAStringSetList(as.list(alleles$alt))
    if (is.null(qual.col)) {
        qual <- rep(NA_real_, nsnp(x))
    } else {
        qual <- getSnpVariable(x, qual.col)
    }
    if (is.null(filter.cols)) {
        filter <- rep(NA_character_, nsnp(x))
    } else {
        filter <- .makeFilter(x, filter.cols)
    }
    S4Vectors::DataFrame(REF=ref, ALT=alt, QUAL=qual, FILTER=filter)
}

.makeInfo <- function(x, info.cols) {
    if (is.null(info.cols)) {
        return(S4Vectors:::make_zero_col_DataFrame(nsnp(x)))
    }
    info.df <- S4Vectors::DataFrame(getSnpVariable(x, info.cols))
    names(info.df) <- info.cols

    type.map <- c(integer="Integer", numeric="Float", character="String",
                  factor="String", logical="Flag")
    meta <- do.call(rbind,
                    lapply(info.cols, function(c) {
                        S4Vectors::DataFrame(Number=if(is.logical(info.df[[c]])) 0 else 1L,
                                             Type=unname(type.map[class(info.df[[c]])]),
                                             Description=.getMetadata(x, c),
                                             row.names=c)
                    }))

    S4Vectors::metadata(info.df)$header <- meta
    info.df
}

.makeGeno <- function(x, ref.allele, scan.index) {
    geno <- getGenotype(x, snp=c(1,-1), scan=c(1,-1))
    geno[ref.allele == "B",] <- 2 - geno[ref.allele == "B",]
    geno <- geno[, scan.index, drop=FALSE]
    geno[is.na(geno)] <- "./."
    geno[geno == 2] <- "0/0"
    geno[geno == 1] <- "0/1"
    geno[geno == 0] <- "1/1"
    geno
}

getScanIndex <- function(x, scan.order, scan.exclude) {
    scanID <- getScanID(x)
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
    match(scan.order, scanID)
}

#setMethod("asVCF", "GenotypeData",
genoDataAsVCF <- function(genoData, sample.col="scanID",
                   id.col="snpID", qual.col=NULL, filter.cols=NULL,
                   info.cols=NULL, scan.exclude=NULL, snp.exclude=NULL,
                   scan.order=NULL, ref.allele=NULL)
{
    if (!(requireNamespace("VariantAnnotation"))) {
        stop("please install VariantAnnotation")
    }

    rowRanges <- .makeRowRanges(genoData, id.col, ref.allele)
    fixed <- .makeFixed(genoData, ref.allele, qual.col, filter.cols)
    info <- .makeInfo(genoData, info.cols)

    samples <- as.character(getScanVariable(genoData, sample.col))
    scan.index <- getScanIndex(genoData, scan.order, scan.exclude)
    colData <- S4Vectors::DataFrame(Samples=seq_along(scan.index),
                                   row.names=samples[scan.index])
    
    geno <- .makeGeno(genoData, ref.allele, scan.index)
    
    format.meta <- S4Vectors::DataFrame(Number=1L, Type="String",
                                        Description="Genotype", row.names="GT")
    info.meta <- S4Vectors::metadata(info)$header
    filt.meta <- attr(fixed$FILTER, "header")
    header <- VariantAnnotation::VCFHeader(samples=samples[scan.index])
    VariantAnnotation::geno(header) <- format.meta
    if (!is.null(info.meta)) VariantAnnotation::info(header) <- info.meta
    if (!is.null(filt.meta)) VariantAnnotation::fixed(header) <- IRanges::DataFrameList(FILTER=filt.meta)
    
    vcf <- VariantAnnotation::VCF(rowRanges=rowRanges, colData=colData,
                                  exptData=list(header=header),
                                  fixed=fixed, info=info, geno=S4Vectors::SimpleList(GT=geno))
    if (!is.null(snp.exclude)) {
        vcf <- vcf[!(getSnpID(genoData) %in% snp.exclude)]
    }
    vcf
}
#)
