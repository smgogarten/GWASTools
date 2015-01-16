# Methods for NcdfGenotypeReader

NcdfGenotypeReader <- function(filename, ...) {
  if (missing(filename)) stop("filename is required")
  new("NcdfGenotypeReader", NcdfReader(filename), ...)
}

setValidity("NcdfGenotypeReader",
    function (object) {
      # check that dimensions and variables are as expected
      
      # check that variables are snpID, chromosome, position, scanID, genotype
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
      if (!hasVariable(object, object@genotypeVar)) {
        return(paste("variable", object@genotypeVar, "not found in", object@filename))
      }
      
      # check that chromosome and position have dimension [snp]
      chrDim <- getDimensionNames(object, object@chromosomeVar)
      if (length(chrDim) != 1 | chrDim != object@snpDim) {
        return(paste("variable", object@chromosomeVar, "has incorrect dimension"))
      }
      posDim <- getDimensionNames(object, object@positionVar)
      if (length(posDim) != 1 | posDim != object@snpDim) {
        return(paste("variable", object@positionVar, "has incorrect dimension"))
      }
      # check that scanID has dimension [scan]
      scanDim <- getDimensionNames(object, object@scanIDvar)
      if (length(scanDim) != 1 | scanDim != object@scanDim) {
        return(paste("variable", object@scanIDvar, "has incorrect dimension"))
      }
      # we don't check snpID because it is a coordinate variable
      
      # check that genotype has dimensions [snp,scan]
      genoDim <- getDimensionNames(object, object@genotypeVar)
      if (length(genoDim) != 2 |
          !all(genoDim == c(object@snpDim, object@scanDim))) {
        return(paste("variable", object@genotypeVar, "has incorrect dimensions"))
      }      

      TRUE
    })

           
# snp and scan are vectors of the format c(start, count)
# count = -1 means read entire dimension
# TODO: modify this function to accept indices or logical vectors
setMethod("getVariable",
          signature(object="NcdfGenotypeReader"),
          function(object, varname, ...) {
              callNextMethod(object, varname, ...)
          })


# accessor methods
# index is logical or integer vector of indices to return
setMethod("getSnpID",
          signature(object="NcdfGenotypeReader"),
          function(object, index) {
            var <- getVariable(object, object@snpIDvar)
            if (missing(index)) var else var[index]
          })

# char=TRUE to return character code
setMethod("getChromosome",
          signature(object="NcdfGenotypeReader"),
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
          signature(object="NcdfGenotypeReader"),
          function(object, index) {
            var <- getVariable(object, object@positionVar)
            if (missing(index)) var else var[index]
          })
                        
setMethod("getScanID",
          signature(object="NcdfGenotypeReader"),
          function(object, index) {
            var <- getVariable(object, object@scanIDvar)
            if (missing(index)) var else var[index]
          })
   
## add names to genotype matrix
.addNamesNcdf <- function(object, var, snp, scan) {
    if (is.matrix(var)) {
        dimnames(var) <- list(getSnpID(object, index=snp), getScanID(object, index=scan))
    }
    var
}

setMethod("getGenotype",
          signature(object="NcdfGenotypeReader"),
          function(object, snp=c(1,-1), scan=c(1,-1), use.names=FALSE, ...) {
              var <- getVariable(object, object@genotypeVar,
                                 start=c(snp[1], scan[1]),
                                 count=c(snp[2], scan[2]), ...)
              if (use.names) {
                  snp.ind <- .startCountToIndex(snp[1], snp[2], nsnp(object))
                  scan.ind <- .startCountToIndex(scan[1], scan[2], nscan(object))
                  var <- .addNamesNcdf(object, var, snp=snp.ind, scan=scan.ind)
              }
              var
          })

setMethod("nsnp", "NcdfGenotypeReader",
          function(object) {
            object@handler$dim[[object@snpDim]]$len
          })

setMethod("nscan", "NcdfGenotypeReader",
          function(object) {
            object@handler$dim[[object@scanDim]]$len
          })

setMethod("autosomeCode", "NcdfGenotypeReader",
          function(object) {
            object@autosomeCode
          })
 
setMethod("XchromCode", "NcdfGenotypeReader",
          function(object) {
            object@XchromCode
          })
          
setMethod("YchromCode", "NcdfGenotypeReader",
          function(object) {
            object@YchromCode
          })
              
setMethod("XYchromCode", "NcdfGenotypeReader",
          function(object) {
            object@XYchromCode
          })
              
setMethod("MchromCode", "NcdfGenotypeReader",
          function(object) {
            object@MchromCode
          })
          
