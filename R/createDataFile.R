

.createNcdf <- function(snp.annotation, filename, variables, n.samples, precision,
                        array.name, genome.build) {

    ## Create the netCDF file and load snp annotation
    ## Define dimensions
    snpdim <- dim.def.ncdf("snp","count", snp.annotation$snpID)
    sampledim <- dim.def.ncdf("sample","count",1:n.samples, unlim=TRUE)

    ## Define variables
    varlist <- list("sampleID"=var.def.ncdf("sampleID", "id", dim=sampledim, missval=0, prec="integer"),
                    "position"=var.def.ncdf("position", "bases", dim=snpdim, missval=-1, prec="integer"),
                    "chromosome"=var.def.ncdf("chromosome", "id", dim=snpdim, missval=-1, prec="integer"))
    if ("genotype" %in% variables) {
        varlist[["genotype"]] <- var.def.ncdf("genotype", "num_A_alleles", dim=list(snpdim,sampledim), missval=-1, prec="byte")
    }
    for (v in setdiff(variables, "genotype")) {
        units <- ifelse(v == "quality", "score", "intensity")
        varlist[[v]] <- var.def.ncdf(v, units, dim=list(snpdim, sampledim), missval=-9999, prec=precision)
    }

    ## Create the netCDF file
    genofile <- create.ncdf(filename, varlist)

    ## Add position data
    put.var.ncdf(genofile, varlist[["position"]], snp.annotation$position)
    put.var.ncdf(genofile, varlist[["chromosome"]], snp.annotation$chromosome)

    ## Add attributes
    if (!is.null(array.name)) att.put.ncdf( genofile, 0, "array_name", array.name )
    if (!is.null(genome.build)) att.put.ncdf( genofile, 0, "genome_build", genome.build )

    close.ncdf(genofile)
}

.createGds <- function(snp.annotation, filename, variables, precision,
                       compress) {

    ## define precision for gds
    precision <- ifelse(precision == "double", "float64", "float32")

    ## define dimensions
    n.snps <- nrow(snp.annotation)

    ## create gds file
    gds <- createfn.gds(filename)

    ## add standard variables
    add.gdsn(gds, "sample.id", storage="integer", valdim=0, compress="") # use valdim=0 then append
    add.gdsn(gds, "snp.id", snp.annotation$snpID, compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.chromosome", snp.annotation$chromosome, storage="uint8",
             compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.position", snp.annotation$position, compress=compress, closezip=TRUE)
    if ("rsID" %in% names(snp.annotation)) {
        add.gdsn(gds, "snp.rs.id", snp.annotation$rsID, compress=compress, closezip=TRUE)
    }
    sync.gds(gds)

    ## add selected variables
    if ("genotype" %in% variables) {
        if (all(c("alleleA", "alleleB") %in% names(snp.annotation))) {
            add.gdsn(gds, "snp.allele", paste(snp.annotation$alleleA, snp.annotation$alleleB, sep="/"),
                     compress=compress, closezip=TRUE)
        }
        geno.node <- add.gdsn(gds, "genotype", storage="bit2", valdim=c(n.snps, 0))
        put.attr.gdsn(geno.node, "snp.order")
    }
    for (v in setdiff(variables, "genotype")) {
        add.gdsn(gds, v, storage=precision, valdim=c(n.snps, 0), compress="")
    }

    sync.gds(gds)
    closefn.gds(gds)
    cleanup.gds(filename)
}

createDataFile <- function(snp.annotation,
                           filename,
                           file.type = c("gds", "ncdf"),
                           variables = "genotype",
                           n.samples = 10,
                           precision = "single",
                           compress = "ZIP.max",
                           array.name = NULL,
                           genome.build = NULL) {

    ## get file type
    file.type <- match.arg(file.type)

    ## check variables
    stopifnot(all(variables %in% c("genotype", "quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio")))

    ## check that snp annotation has right columns
    stopifnot(all(c("snpID", "chromosome", "position") %in% names(snp.annotation)))

    ## make sure all snp annotation columns are integers
    if (!is(snp.annotation$snpID, "integer")) {
        snp.annotation$snpID <- as.integer(snp.annotation$snpID)
        warning(paste("coerced snpID to type integer"))
    }
    if (!is(snp.annotation$chromosome, "integer")) {
        snp.annotation$chromosome <- as.integer(snp.annotation$chromosome)
        warning(paste("coerced chromosome to type integer"))
    }
    if (!is(snp.annotation$position, "integer")) {
        snp.annotation$position <- as.integer(snp.annotation$position)
        warning(paste("coerced position to type integer"))
    }

    ## make sure snpID is unique
    stopifnot(length(snp.annotation$snpID) == length(unique(snp.annotation$snpID)))

    ## make sure snpID is sorted by chromsome and position
    stopifnot(all(snp.annotation$snpID == sort(snp.annotation$snpID)))
    sorted <- order(snp.annotation$chromosome, snp.annotation$position)
    if (!all(snp.annotation$snpID == snp.annotation$snpID[sorted])) {
        stop("snpID is not sorted by chromosome and position")
    }

    if (file.type == "gds") {
        ## don't need n.samples since we will use append later
        .createGds(snp.annotation, filename, variables, precision, compress)
    } else if (file.type == "ncdf") {
        .createNcdf(snp.annotation, filename, variables, n.samples, precision, array.name, genome.build)
    }

    return(invisible(NULL))
}

