# Methods for GdsIntensityReader

GdsIntensityReader <- function(filename, ...) {
  if (missing(filename)) stop("filename is required")
  
  input.gds <- is(filename, 'gds.class')
  tmpobj <- GdsReader(filename)
  
  tryCatch(new("GdsIntensityReader", tmpobj, ...),
           error=function(e) {if (!input.gds) close(tmpobj)
             stop(e)})
}

setValidity("GdsIntensityReader",
    function (object) {
      # check that dimensions and variables are as expected

      # check that variables are snpID, chromosome, position, scanID
      if (!hasVariable(object, object@snpIDvar)) {
        return(paste("variable", object@snpIDvar, "not found in", object@filename))
      }
      if (!hasVariable(object, object@chromosomeVar)) {
        return(paste("variable", object@chromosomeVar, "not found in", object@filename))
      }
      if (!hasVariable(object, object@positionVar)) {
        return(paste("variable", object@positionVar, "not found in", object@filename))
      }
      if (!hasVariable(object, object@scanIDvar)) {
        return(paste("variable", object@scanIDvar, "not found in", object@filename))
      }

      # check that chromosome and position have same dimension as snpID
      snpDim <- getDimension(object, object@snpIDvar)
      chrDim <- getDimension(object, object@chromosomeVar)
      if (chrDim != snpDim) {
        return(paste("variable", object@chromosomeVar, "has incorrect dimension"))
      }
      posDim <- getDimension(object, object@positionVar)
      if (posDim != snpDim) {
        return(paste("variable", object@positionVar, "has incorrect dimension"))
      }

      # check that other variables have dimensions [snp,scan]
      scanDim <- getDimension(object, object@scanIDvar)
      for (var in c(object@qualityVar, object@xVar, object@yVar,
                    object@bafVar, object@lrrVar)) {
        if (hasVariable(object, var)) {
            varDim <- getDimension(object, var)
            if (length(varDim) != 2 |
                !all(varDim == c(snpDim, scanDim))) {
                return(paste("variable", var, "has incorrect dimensions"))
            }
        }
      }

      TRUE
    })


# snp and scan are vectors of the format c(start, count)
# count = -1 means read entire dimension
# TODO: modify this function to accept indices or logical vectors
setMethod("getVariable",
          signature(object="GdsIntensityReader"),
          function(object, varname, snp, scan, ...) {

            if (!missing(snp) & !missing(scan)) {
              # get start and count from snp
              snpstart = snp[1]
              snpcount = snp[2]

              # get start and count from scan
              scanstart = scan[1]
              scancount = scan[2]

              callNextMethod(object, varname, start=c(snpstart, scanstart),
                             count=c(snpcount, scancount), ...)
            } else {
              callNextMethod(object, varname, ...)
            }
          })


# accessor methods
# index is logical or integer vector of indices to return
setMethod("getSnpID",
          signature(object="GdsIntensityReader"),
          function(object, index) {
            var <- getVariable(object, object@snpIDvar)
            if (missing(index)) var else var[index]
          })

# char=TRUE to return character code
setMethod("getChromosome",
          signature(object="GdsIntensityReader"),
          function(object, index, char=FALSE) {
            var <- getVariable(object, object@chromosomeVar)
            if (!missing(index)) var <- var[index]

            # convert to characters
            if (char) {
              # default is unknown code
              chromChar <- rep("U", length(var))
              autosome <- var %in% object@autosomeCode
              chromChar[autosome] <- as.character(var[autosome])
              xchrom <- var == object@XchromCode & !is.na(var)
              chromChar[xchrom] <- "X"
              ychrom <- var == object@YchromCode & !is.na(var)
              chromChar[ychrom] <- "Y"
              xychrom <- var == object@XYchromCode & !is.na(var)
              chromChar[xychrom] <- "XY"
              mchrom <- var == object@MchromCode & !is.na(var)
              chromChar[mchrom] <- "M"
              var <- chromChar
            }
            var
          })

setMethod("getPosition",
          signature(object="GdsIntensityReader"),
          function(object, index) {
            var <- getVariable(object, object@positionVar)
            if (missing(index)) var else var[index]
          })

setMethod("getScanID",
          signature(object="GdsIntensityReader"),
          function(object, index) {
            var <- getVariable(object, object@scanIDvar)
            if (missing(index)) var else var[index]
          })

setMethod("hasQuality",
          signature(object="GdsIntensityReader"),
          function(object) {
            hasVariable(object, object@qualityVar)
          })

setMethod("getQuality",
          signature(object="GdsIntensityReader"),
          function(object, ...) {
            getVariable(object, object@qualityVar, ...)
          })

setMethod("hasX",
          signature(object="GdsIntensityReader"),
          function(object) {
            hasVariable(object, object@xVar)
          })

setMethod("getX",
          signature(object="GdsIntensityReader"),
          function(object, ...) {
            getVariable(object, object@xVar, ...)
          })

setMethod("hasY",
          signature(object="GdsIntensityReader"),
          function(object) {
            hasVariable(object, object@yVar)
          })

setMethod("getY",
          signature(object="GdsIntensityReader"),
          function(object, ...) {
            getVariable(object, object@yVar, ...)
          })

setMethod("hasBAlleleFreq",
          signature(object="GdsIntensityReader"),
          function(object) {
            hasVariable(object, object@bafVar)
          })

setMethod("getBAlleleFreq",
          signature(object="GdsIntensityReader"),
          function(object, ...) {
            getVariable(object, object@bafVar, ...)
          })

setMethod("hasLogRRatio",
          signature(object="GdsIntensityReader"),
          function(object) {
            hasVariable(object, object@lrrVar)
          })

setMethod("getLogRRatio",
          signature(object="GdsIntensityReader"),
          function(object, ...) {
            getVariable(object, object@lrrVar, ...)
          })

setMethod("nsnp", "GdsIntensityReader",
          function(object) {
            getDimension(object, object@snpIDvar)
          })

setMethod("nscan", "GdsIntensityReader",
          function(object) {
            getDimension(object, object@scanIDvar)
          })

setMethod("autosomeCode", "GdsIntensityReader",
          function(object) {
            object@autosomeCode
          })

setMethod("XchromCode", "GdsIntensityReader",
          function(object) {
            object@XchromCode
          })

setMethod("YchromCode", "GdsIntensityReader",
          function(object) {
            object@YchromCode
          })

setMethod("XYchromCode", "GdsIntensityReader",
          function(object) {
            object@XYchromCode
          })

setMethod("MchromCode", "GdsIntensityReader",
          function(object) {
            object@MchromCode
          })

