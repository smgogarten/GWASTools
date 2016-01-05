
.createNcdf <- function(snp.annotation, filename, variables, n.samples,
                        precision="single",
                        array.name=NULL, genome.build=NULL,
                        var.data=NULL) {

    ## Create the netCDF file and load snp annotation
    ## Define dimensions
    snpdim <- ncdim_def("snp","count", snp.annotation$snpID)
    sampledim <- ncdim_def("sample","count",1:n.samples, unlim=TRUE)

    ## Define variables
    varlist <- list("sampleID"=ncvar_def("sampleID", "id", dim=sampledim, missval=0, prec="integer"),
                    "position"=ncvar_def("position", "bases", dim=snpdim, missval=-1, prec="integer"),
                    "chromosome"=ncvar_def("chromosome", "id", dim=snpdim, missval=-1, prec="integer"))
    if ("genotype" %in% variables) {
        varlist[["genotype"]] <- ncvar_def("genotype", "num_A_alleles", dim=list(snpdim,sampledim), missval=-1, prec="byte")
    }
    for (v in setdiff(variables, "genotype")) {
        units <- ifelse(v == "quality", "score", "intensity")
        varlist[[v]] <- ncvar_def(v, units, dim=list(snpdim, sampledim), missval=-9999, prec=precision)
    }

    ## Create the netCDF file
    genofile <- nc_create(filename, varlist)

    ## Add position data
    ncvar_put(genofile, varlist[["position"]], snp.annotation$position)
    ncvar_put(genofile, varlist[["chromosome"]], snp.annotation$chromosome)

    ## Add any other variable data
    if (!is.null(var.data)) {
        stopifnot(all(names(var.data) %in% c("sampleID", variables)))
        for (v in names(var.data)) {
            ncvar_put(genofile, varlist[[v]], var.data[[v]])
        }
    }

    ## Add attributes
    if (!is.null(array.name)) ncatt_put( genofile, 0, "array_name", array.name )
    if (!is.null(genome.build)) ncatt_put( genofile, 0, "genome_build", genome.build )

    nc_sync(genofile)
    genofile
}

.createGds <- function(snp.annotation, filename, variables, precision="single",
                       compress="ZIP_RA:8M", compress.geno="", compress.annot="ZIP_RA",
                       sample.storage="integer", var.data=NULL) {

    ## define precision for gds
    precision <- ifelse(precision == "double", "float64", "float32")

    ## define dimensions
    n.snps <- nrow(snp.annotation)

    ## create gds file
    gds <- createfn.gds(filename)

    ## add standard variables
    add.gdsn(gds, "sample.id", storage=sample.storage, valdim=0, compress=compress.annot) # use valdim=0 then append
    add.gdsn(gds, "snp.id", snp.annotation$snpID, compress=compress.annot, closezip=TRUE)
    add.gdsn(gds, "snp.chromosome", snp.annotation$chromosome, storage="uint8",
             compress=compress.annot, closezip=TRUE)
    add.gdsn(gds, "snp.position", snp.annotation$position, compress=compress.annot, closezip=TRUE)
    if ("snpName" %in% names(snp.annotation))
        add.gdsn(gds, "snp.rs.id", snp.annotation$snpName, compress=compress.annot, closezip=TRUE)
    sync.gds(gds)

    ## add selected variables
    if ("genotype" %in% variables) {
        if (all(c("alleleA", "alleleB") %in% names(snp.annotation))) {
            add.gdsn(gds, "snp.allele", paste(snp.annotation$alleleA, snp.annotation$alleleB, sep="/"),
                     compress=compress.annot, closezip=TRUE)
        }
        geno.node <- add.gdsn(gds, "genotype", storage="bit2", valdim=c(n.snps, 0), compress=compress.geno)
        put.attr.gdsn(geno.node, "snp.order")
    }
    for (v in setdiff(variables, "genotype")) {
        add.gdsn(gds, v, storage=precision, valdim=c(n.snps, 0), compress=compress)
    }

    ## Add any other variable data
    if (!is.null(var.data)) {
        stopifnot(all(names(var.data) %in% c("sample.id", variables)))
        for (v in names(var.data)) {
            append.gdsn(index.gdsn(gds, v), var.data[[v]])
        }
    }

    sync.gds(gds)
    gds
}

.createGdsBySnp <- function(sample.id, snp.annotation, filename, variables, precision,
                       compress) {

    ## define precision for gds
    precision <- ifelse(precision == "double", "float64", "float32")

    ## create gds file
    gds <- createfn.gds(filename)

    ## add standard variables
    add.gdsn(gds, "sample.id", sample.id, compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.id", snp.annotation$snpID, compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.chromosome", snp.annotation$chromosome, storage="uint8",
             compress=compress, closezip=TRUE)
    add.gdsn(gds, "snp.position", snp.annotation$position, compress=compress, closezip=TRUE)
    sync.gds(gds)

    ## add selected variables
    for (v in variables) {
        add.gdsn(gds, v, storage=precision, valdim=c(nrow(snp.annotation), length(sample.id)))
    }

    sync.gds(gds)
    gds
}

.addData <- function(x, ...) UseMethod(".addData", x)
.addData.gds.class <- function(x, vars, dat, sample.id, ...) {
    append.gdsn(index.gdsn(x, "sample.id"), val=sample.id)
    for (v in vars) {
        ## set missing code for genotype
        if (v == "genotype") dat[[v]][is.na(dat[[v]])] <- 3
        append.gdsn(index.gdsn(x, v), val=dat[[v]])
    }
}

.addData.ncdf4 <- function(x, vars, dat, sample.id, sample.index) {
    ncvar_put(x, "sampleID", vals=sample.id, start=sample.index, count=1)
    for (v in vars) {
        ## set missing code for genotype
        if (v == "genotype") dat[[v]][is.na(dat[[v]])] <- -1
        ncvar_put(x, v, vals=dat[[v]], start=c(1,sample.index), count=c(-1,1))
    }
}

.addDataBySnp <- function(x, ...) UseMethod(".addDataBySnp", x)
.addDataBySnp.gds.class <- function(x, vars, dat, snp.start, snp.count) {
    for (v in vars) {
        write.gdsn(index.gdsn(x, v), val=dat[[v]], start=c(snp.start,1), count=c(snp.count,-1))
    }
}

.addDataBySnp.ncdf4 <- function(x, vars, dat, snp.start, snp.count) {
    for (v in vars) {
        ncvar_put(x, v, vals=dat[[v]], start=c(snp.start,1), count=c(snp.count,-1))
    }
}

.close <- function(x, ...) UseMethod(".close", x)
.close.gds.class <- function(x, verbose=FALSE) {
    vars <- ls.gdsn(x)
    vars <- vars[!grepl("^snp", vars)] # snp nodes already done
    for (v in vars) readmode.gdsn(index.gdsn(x, v))
    sync.gds(x)

    ## close and cleanup
    filename <- x$filename
    closefn.gds(x)
    cleanup.gds(filename, verbose=verbose)
}
.close.ncdf4 <- function(x, ...) nc_close(x)



.createGdsDosage <- function(snp.df, scan.df, filename, genotypeDim, miss.val, precision="single",
                             compress.annot="ZIP_RA") {

  # redefine precision for gds
  precision <- ifelse(precision == "double", "float64",
                      ifelse(precision == "single", "float32", precision))

  # create GDS
  gfile <- createfn.gds(filename)

  add.gdsn(gfile, "snp.id", snp.df$snpID, compress=compress.annot, closezip=TRUE)
  add.gdsn(gfile, "sample.id", scan.df$scanID, compress=compress.annot, closezip=TRUE)

  n <- add.gdsn(gfile, name="description")
  put.attr.gdsn(n, "FileFormat", "IMPUTED_DOSAGE")

  geno.valdim <- switch(genotypeDim,
                        "snp,scan"=c(length(snp.df$snpID), length(scan.df$scanID)),
                        "scan,snp"=c(length(scan.df$scanID), length(snp.df$snpID)))
  gGeno <- add.gdsn(gfile, "genotype", valdim=geno.valdim, storage=precision)

  geno.order <- switch(genotypeDim,
                       "snp,scan"="snp.order",
                       "scan,snp"="sample.order")
  put.attr.gdsn(gGeno, geno.order)
  put.attr.gdsn(gGeno, "missing.value", miss.val)

  sync.gds(gfile)
  gfile
}

.createNcdfDosage <- function(snp.df, scan.df, filename, miss.val, precision="single") {
  # define dimensions
  snpdim <- ncdim_def("snp", "count", snp.df$snpID)
  sampledim <- ncdim_def("sample", "count", scan.df$scanID, unlim=TRUE)
  chardim <- ncdim_def("nchar", "", 1)

  # define variables
  varID <- ncvar_def("sampleID", "id", dim=sampledim, missval=0, prec="integer")
  varpos <- ncvar_def("position", "bases", dim=snpdim, missval=-1, prec="integer")
  varchr <- ncvar_def("chromosome", "id", dim=snpdim, missval=-1, prec="integer")
  varA <- ncvar_def("alleleA", "allele", dim=list(chardim,snpdim), missval="0", prec="char")
  varB <- ncvar_def("alleleB", "allele", dim=list(chardim,snpdim), missval="0", prec="char")
  vargeno <- ncvar_def("genotype", "A_allele_dosage", dim=list(snpdim,sampledim), missval=miss.val, prec=precision)

  # create the NetCDF file
  nc <- nc_create(filename, list(varID, varpos, varchr, varA, varB, vargeno))

  ncvar_put(nc, varID, scan.df$scanID)

  nc
}


.addDosage <- function(x, ...) UseMethod(".addDosage", x)
.addDosage.gds.class <- function(x, dosage, start, count) {
    write.gdsn(index.gdsn(x, "genotype"), dosage, start=start, count=count)
}

.addDosage.ncdf4 <- function(x, dosage, start, count) {
    ncvar_put(x, "genotype", dosage, start=start, count=count)
}

.addSnpVars <- function(x, ...) UseMethod(".addSnpVars", x)
.addSnpVars.gds.class <- function(x, snpAnnot, compress) {
  add.gdsn(x, "snp.chromosome", snpAnnot$chromosome, compress=compress, closezip=TRUE)
  add.gdsn(x, "snp.position", snpAnnot$position, compress=compress, closezip=TRUE)
  add.gdsn(x, "snp.allele", paste(snpAnnot$alleleA, snpAnnot$alleleB, sep="/"), compress=compress, closezip=TRUE)

}
.addSnpVars.ncdf4 <- function(x, snpAnnot, ...) {
  ncvar_put(x, "position", snpAnnot$position)
  ncvar_put(x, "chromosome", snpAnnot$chromosome)
  ncvar_put(x, "alleleA", snpAnnot$alleleA)
  ncvar_put(x, "alleleB", snpAnnot$alleleB)
}

.reopenGds <- function(x) {
    closefn.gds(x)
    openfn.gds(x$filename, readonly=FALSE)
}

.compressDosage <- function(x, compress) {
    compression.gdsn(index.gdsn(x, "genotype"), compress=compress)
}

.addChromosomeAttributes <- function(new, old) {
    chr.node <- index.gdsn(new, "snp.chromosome")
    put.attr.gdsn(chr.node, "autosome.start", min(autosomeCode(old)))
    put.attr.gdsn(chr.node, "autosome.end", max(autosomeCode(old)))
    put.attr.gdsn(chr.node, "X", XchromCode(old))
    put.attr.gdsn(chr.node, "XY", XYchromCode(old))
    put.attr.gdsn(chr.node, "Y", YchromCode(old))
    put.attr.gdsn(chr.node, "M", MchromCode(old))
    put.attr.gdsn(chr.node, "MT", MchromCode(old))
}
