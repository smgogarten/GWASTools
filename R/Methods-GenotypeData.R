# Methods for GenotypeData

# constructor
GenotypeData <- function(data, snpAnnot=NULL, scanAnnot=NULL) {
  new("GenotypeData", data=data, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
}

# validity method
setValidity("GenotypeData",
          function(object) {
            # if snpAnnot is given, check
            if (!is.null(object@snpAnnot)) {

              # check that codes match
              if (XchromCode(object@snpAnnot) != XchromCode(object@data) |
                  YchromCode(object@snpAnnot) != YchromCode(object@data) |
                  XYchromCode(object@snpAnnot) != XYchromCode(object@data) |
                  MchromCode(object@snpAnnot) != MchromCode(object@data)) {
                return("data and snpAnnot have different chromosome codes")
              }
              
              # check that snpIDs match
              snpID <- getSnpID(object@snpAnnot)
              dataID <- getSnpID(object@data)

              # check length
              if (length(snpID) != length(dataID)) {
                return("data and snpAnnot have different lengths")
              }

              # check elements
              if (!allequal(snpID, dataID)) {
                return("data and snpAnnot have different snpIDs")
              }

              # check that chromosomes match
              snpChr <- getChromosome(object@snpAnnot)
              dataChr <- getChromosome(object@data)
              if (!allequal(snpChr, dataChr)) {
                return("data and snpAnnot have different chromosomes")
              }
              
              # check that positions match
              snpPos <- getPosition(object@snpAnnot)
              dataPos <- getPosition(object@data)
              if (!allequal(snpPos, dataPos)) {
                return("data and snpAnnot have different positions")
              }
            }

            #if scanAnnot is given, check
            if (!is.null(object@scanAnnot)) {
              
              # check that scanIDs match
              scanID <- getScanID(object@scanAnnot)
              dataID <- getScanID(object@data)

              # check length
              if (length(scanID) != length(dataID)) {
                return("data and scanAnnot have different lengths")
              }

              # check elements
              if (!allequal(scanID, dataID)) {
                return("data and scanAnnot have different scanIDs")
              }
            }
            TRUE
          })

setMethod("show",
          signature(object = "GenotypeData"),
          function(object) {
            cat("An object of class", class(object), "\n")
            cat(" | data:\n")
            show(object@data)
            cat(" | SNP Annotation:\n")
            show(object@snpAnnot)
            cat(" | Scan Annotation:\n")
            show(object@scanAnnot)
          })

setMethod("hasSnpAnnotation",
          signature(object = "GenotypeData"),
          function(object) {
            !is.null(object@snpAnnot)
          })

setMethod("hasScanAnnotation",
          signature(object = "GenotypeData"),
          function(object) {
            !is.null(object@scanAnnot)
          })

setMethod("getSnpID",
          signature(object = "GenotypeData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getSnpID(object@snpAnnot, ...)
            } else {
              getSnpID(object@data, ...)
            }
          })

setMethod("getChromosome",
          signature(object = "GenotypeData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getChromosome(object@snpAnnot, ...)
            } else {
              getChromosome(object@data, ...)
            }
          })

setMethod("getPosition",
          signature(object = "GenotypeData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getPosition(object@snpAnnot, ...)
            } else {
              getPosition(object@data, ...)
            }
          })

setMethod("getAlleleA",
          signature(object = "GenotypeData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getAlleleA(object@snpAnnot, ...)
            } else {
              NULL
            }
          })

setMethod("getAlleleB",
          signature(object = "GenotypeData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getAlleleB(object@snpAnnot, ...)
            } else {
              NULL
            }
          })

setMethod("getScanID",
          signature(object = "GenotypeData"),
          function(object, ...) {
            if (hasScanAnnotation(object)) {
              getScanID(object@scanAnnot, ...)
            } else {
              getScanID(object@data, ...)
            }
          })

setMethod("hasSex",
          signature(object = "GenotypeData"),
          function(object) {
            if (hasScanAnnotation(object)) {
              hasSex(object@scanAnnot)
            } else FALSE
          })

setMethod("getSex",
          signature(object = "GenotypeData"),
          function(object, ...) {
            if (hasScanAnnotation(object)) {
              getSex(object@scanAnnot, ...)
            } else {
              warning("scan annotation not found")
              return(NULL)
            }
          })

setMethod("hasSnpVariable",
          signature(object = "GenotypeData"),
          function(object, varname) {
            if (hasSnpAnnotation(object)) {
              hasVariable(object@snpAnnot, varname)
            } else FALSE
          })

setMethod("getSnpVariable",
          signature(object = "GenotypeData"),
          function(object, varname, ...) {
            if (hasSnpAnnotation(object)) {
              getVariable(object@snpAnnot, varname, ...)
            } else {
              warning("snp annotation not found")
              return(NULL)
            }
          })

setMethod("getSnpVariableNames",
          signature(object = "GenotypeData"),
          function(object) {
            if (hasSnpAnnotation(object)) {
              getVariableNames(object@snpAnnot)
            } else {
              warning("snp annotation not found")
              return(NULL)
            }
          })

setMethod("hasScanVariable",
          signature(object = "GenotypeData"),
          function(object, varname) {
            if (hasScanAnnotation(object)) {
              hasVariable(object@scanAnnot, varname)
            } else FALSE
          })

setMethod("getScanVariable",
          signature(object = "GenotypeData"),
          function(object, varname, ...) {
            if (hasScanAnnotation(object)) {
              getVariable(object@scanAnnot, varname, ...)
            } else {
              warning("scan annotation not found")
              return(NULL)
            }
          })

setMethod("getScanVariableNames",
          signature(object = "GenotypeData"),
          function(object) {
            if (hasScanAnnotation(object)) {
              getVariableNames(object@scanAnnot)
            } else {
              warning("scan annotation not found")
              return(NULL)
            }
          })

setMethod("hasVariable",
          signature(object = "GenotypeData"),
          function(object, varname) {
            hasVariable(object@data, varname)
          })

setMethod("getVariable",
          signature(object = "GenotypeData"),
          function(object, varname, ...) {
            getVariable(object@data, varname, ...)
          })

# char=TRUE to return character values "A/B"
setMethod("getGenotype",
          signature(object = "GenotypeData"),
          function(object, char=FALSE, sort=TRUE, ...) {
            geno <- getGenotype(object@data, ...)
            if (char) {
              alleleA <- getAlleleA(object)
              alleleB <- getAlleleB(object)
              geno <- genotypeToCharacter(geno, alleleA, alleleB, sort)
            }
            geno
          })

setMethod("nsnp", "GenotypeData",
          function(object) {
            nsnp(object@data)
          })

setMethod("nscan", "GenotypeData",
          function(object) {
            nscan(object@data)
          })

setMethod("open",
    signature(con = "GenotypeData"),
    function (con, ...) {
      open(con@data, ...)
    })

setMethod("close",
    signature(con = "GenotypeData"),
    function (con, ...) {
      close(con@data, ...)
    })

setMethod("XchromCode", "GenotypeData",
          function(object) {
            XchromCode(object@data)
          })
          
setMethod("YchromCode", "GenotypeData",
          function(object) {
            YchromCode(object@data)
          })
              
setMethod("XYchromCode", "GenotypeData",
          function(object) {
            XYchromCode(object@data)
          })
              
setMethod("MchromCode", "GenotypeData",
          function(object) {
            MchromCode(object@data)
          })
          
