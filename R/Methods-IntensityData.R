# Methods for IntensityData

# constructor
IntensityData <- function(data, snpAnnot=NULL, scanAnnot=NULL) {
  new("IntensityData", data=data, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
}

# validity method
setValidity("IntensityData",
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
          signature(object = "IntensityData"),
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
          signature(object = "IntensityData"),
          function(object) {
            !is.null(object@snpAnnot)
          })

setMethod("hasScanAnnotation",
          signature(object = "IntensityData"),
          function(object) {
            !is.null(object@scanAnnot)
          })

setMethod("getSnpID",
          signature(object = "IntensityData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getSnpID(object@snpAnnot, ...)
            } else {
              getSnpID(object@data, ...)
            }
          })

setMethod("getChromosome",
          signature(object = "IntensityData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getChromosome(object@snpAnnot, ...)
            } else {
              getChromosome(object@data, ...)
            }
          })

setMethod("getPosition",
          signature(object = "IntensityData"),
          function(object, ...) {
            if (hasSnpAnnotation(object)) {
              getPosition(object@snpAnnot, ...)
            } else {
              getPosition(object@data, ...)
            }
          })

setMethod("getScanID",
          signature(object = "IntensityData"),
          function(object, ...) {
            if (hasScanAnnotation(object)) {
              getScanID(object@scanAnnot, ...)
            } else {
              getScanID(object@data, ...)
            }
          })

setMethod("hasSex",
          signature(object = "IntensityData"),
          function(object) {
            if (hasScanAnnotation(object)) {
              hasSex(object@scanAnnot)
            } else FALSE
          })

setMethod("getSex",
          signature(object = "IntensityData"),
          function(object, ...) {
            if (hasScanAnnotation(object)) {
              getSex(object@scanAnnot, ...)
            } else {
              warning("scan annotation not found")
              return(NULL)
            }
          })

setMethod("hasSnpVariable",
          signature(object = "IntensityData"),
          function(object, varname) {
            if (hasSnpAnnotation(object)) {
              hasVariable(object@snpAnnot, varname)
            } else FALSE
          })

setMethod("getSnpVariable",
          signature(object = "IntensityData"),
          function(object, varname, ...) {
            if (hasSnpAnnotation(object)) {
              getVariable(object@snpAnnot, varname, ...)
            } else {
              warning("snp annotation not found")
              return(NULL)
            }
          })

setMethod("getSnpVariableNames",
          signature(object = "IntensityData"),
          function(object) {
            getVariableNames(object@snpAnnot)
          })

setMethod("hasScanVariable",
          signature(object = "IntensityData"),
          function(object, varname) {
            if (hasScanAnnotation(object)) {
              hasVariable(object@scanAnnot, varname)
            } else FALSE
          })

setMethod("getScanVariable",
          signature(object = "IntensityData"),
          function(object, varname, ...) {
            if (hasScanAnnotation(object)) {
              getVariable(object@scanAnnot, varname, ...)
            } else {
              warning("scan annotation not found")
              return(NULL)
            }
          })

setMethod("getScanVariableNames",
          signature(object = "IntensityData"),
          function(object) {
            getVariableNames(object@scanAnnot)
          })

setMethod("hasVariable",
          signature(object = "IntensityData"),
          function(object, varname) {
            hasVariable(object@data, varname)
          })

setMethod("getVariable",
          signature(object = "IntensityData"),
          function(object, varname, ...) {
            getVariable(object@data, varname, ...)
          })

setMethod("hasQuality",
          signature(object = "IntensityData"),
          function(object) {
            hasQuality(object@data)
          })

setMethod("getQuality",
          signature(object = "IntensityData"),
          function(object, ...) {
            getQuality(object@data, ...)
          })

setMethod("hasX",
          signature(object = "IntensityData"),
          function(object) {
            hasX(object@data)
          })

setMethod("getX",
          signature(object = "IntensityData"),
          function(object, ...) {
            getX(object@data, ...)
          })

setMethod("hasY",
          signature(object = "IntensityData"),
          function(object) {
            hasY(object@data)
          })

setMethod("getY",
          signature(object = "IntensityData"),
          function(object, ...) {
            getY(object@data, ...)
          })

setMethod("hasBAlleleFreq",
          signature(object = "IntensityData"),
          function(object) {
            hasBAlleleFreq(object@data)
          })

setMethod("getBAlleleFreq",
          signature(object = "IntensityData"),
          function(object, ...) {
            getBAlleleFreq(object@data, ...)
          })

setMethod("hasLogRRatio",
          signature(object = "IntensityData"),
          function(object) {
            hasLogRRatio(object@data)
          })

setMethod("getLogRRatio",
          signature(object = "IntensityData"),
          function(object, ...) {
            getLogRRatio(object@data, ...)
          })


setMethod("nsnp", "IntensityData",
          function(object) {
            nsnp(object@data)
          })

setMethod("nscan", "IntensityData",
          function(object) {
            nscan(object@data)
          })

setMethod("open",
    signature(con = "IntensityData"),
    function (con, ...) {
      open(con@data, ...)
    })

setMethod("close",
    signature(con = "IntensityData"),
    function (con, ...) {
      close(con@data, ...)
    })

setMethod("XchromCode", "IntensityData",
          function(object) {
            XchromCode(object@data)
          })
          
setMethod("YchromCode", "IntensityData",
          function(object) {
            YchromCode(object@data)
          })
              
setMethod("XYchromCode", "IntensityData",
          function(object) {
            XYchromCode(object@data)
          })
              
setMethod("MchromCode", "IntensityData",
          function(object) {
            MchromCode(object@data)
          })
          


