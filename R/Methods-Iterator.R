GenotypeIterator <- function(genoData, snpFilter) {
    class(genoData) <- "GenotypeIterator"
    genoData@snpFilter <- snpFilter

    ## pass-by-reference slot for lastFilter
    genoData@lastFilter <- new.env()
    lastFilter(genoData) <- 1

    genoData
}


GenotypeBlockIterator <- function(genoData, snpBlock=10000, snpInclude=NULL) {
    snps <- 1:nsnp(genoData)
    if (!is.null(snpInclude)) {
        snps <- which(getSnpID(genoData) %in% snpInclude)
    }
    
    if (snpBlock > length(snps)) {
        snpBlock <- length(snps)
    }

    snpFilter <- split(snps, ceiling((1:length(snps))/snpBlock))

    genoData <- GenotypeIterator(genoData, snpFilter)
    class(genoData) <- "GenotypeBlockIterator"
    genoData@snpBlock <- as.integer(snpBlock)

    genoData
}


setMethod("snpFilter",
          "GenotypeIterator",
          function(x) {
              x@snpFilter
          })

setMethod("lastFilter",
          "GenotypeIterator",
          function(x) {
              x@lastFilter$i
          })

setReplaceMethod("lastFilter",
                 c("GenotypeIterator", "numeric"),
                 function(x, value) {
                     x@lastFilter$i <- as.integer(value)
                     x
                 })

setMethod("currentFilter",
          "GenotypeIterator",
          function(x) {
              snpFilter(x)[[lastFilter(x)]]
          })

setMethod("iterateFilter",
          "GenotypeIterator",
          function(x) {
              ## set filter for next element
              if (lastFilter(x) < length(snpFilter(x))) {
                  i <- lastFilter(x) + 1
                  lastFilter(x) <- i
                  return(TRUE)
              } else {
                  return(FALSE)
              }
          })

setMethod("getSnpID",
          "GenotypeIterator",
          function(object, ...) {
              callNextMethod(object, index=currentFilter(object), ...)
          })

setMethod("getChromosome",
          "GenotypeIterator",
          function(object, ...) {
              callNextMethod(object, index=currentFilter(object), ...)
          })

setMethod("getPosition",
          "GenotypeIterator",
          function(object, ...) {
              callNextMethod(object, index=currentFilter(object), ...)
          })

setMethod("getAlleleA",
          "GenotypeIterator",
          function(object, ...) {
              callNextMethod(object, index=currentFilter(object), ...)
          })

setMethod("getAlleleB",
          "GenotypeIterator",
          function(object, ...) {
              callNextMethod(object, index=currentFilter(object), ...)
          })

setMethod("getSnpVariable",
          "GenotypeIterator",
          function(object, varname, ...) {
              callNextMethod(object, varname=varname, index=currentFilter(object), ...)
          })

setMethod("getGenotypeSelection",
          "GenotypeIterator",
          function(object, ...) {
              callNextMethod(object, snp=currentFilter(object), ...)
          })
