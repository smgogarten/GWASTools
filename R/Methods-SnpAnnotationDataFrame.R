# Methods for SnpAnnotationDataFrame

# constructor
SnpAnnotationDataFrame <- function(data, metadata, ...) {
  if (missing(metadata)) {
    new("SnpAnnotationDataFrame", data=data, dimLabels=c("snps", "variables"), ...)
  } else {
    new("SnpAnnotationDataFrame", data=data, varMetadata=metadata,
        dimLabels=c("snps", "variables"), ...)
  }
}

# validity
setValidity("SnpAnnotationDataFrame",
          function(object) {
            # check for required columns
            if (!hasVariable(object, object@idCol)) {
              return(paste("missing required column", object@idCol))
            }
            if (!hasVariable(object, object@chromosomeCol)) {
              return(paste("missing required column", object@chromosomeCol))
            }
            if (!hasVariable(object, object@positionCol)) {
              return(paste("missing required column", object@positionCol))
            }

            # check format of required columns
            # snpID should be a unique integer vector
            snpID <- getSnpID(object)
            if (length(snpID) != length(unique(snpID)) ||
                !is.integer(snpID)) {
              return(paste(object@idCol, "must be a unique integer vector"))
            }
            # chromosome should be an integer vector
            if (!is.integer(getChromosome(object))) {
              return(paste(object@chromosomeCol, "must be an integer vector"))
            }
            # position should be an integer vector
            if (!is.integer(getPosition(object))) {
              return(paste(object@positionCol, "must be an integer vector"))
            }
            TRUE
          })
        
setMethod("hasVariable",
          signature(object = "SnpAnnotationDataFrame"),
          function(object, varname) {
            varname %in% varLabels(object)
          })
            
setMethod("getVariable",
          signature(object = "SnpAnnotationDataFrame"),
          function(object, varname, index) {
            # check that variable exists
            varexist <- hasVariable(object, varname)
            if (!all(varexist)) {
              warning(paste(paste(varname[!varexist], collapse=","), "not found"))
              return(NULL)
            }
            
            if (missing(index)) {
              object@data[, varname]
            } else {
              object@data[index, varname]
            }
          })
 
setMethod("getSnpID",
          signature(object = "SnpAnnotationDataFrame"),
          function(object, ...) {
            getVariable(object, object@idCol, ...)
          })
        
# char=TRUE to return character code      
setMethod("getChromosome",
          signature(object = "SnpAnnotationDataFrame"),
          function(object, char=FALSE, ...) {
            var <- getVariable(object, object@chromosomeCol, ...)
            
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
          signature(object = "SnpAnnotationDataFrame"),
          function(object, ...) {
            getVariable(object, object@positionCol, ...)
          })
      
setMethod("getAlleleA",
          signature(object = "SnpAnnotationDataFrame"),
          function(object, ...) {
            getVariable(object, object@alleleACol, ...)
          })

setMethod("getAlleleB",
          signature(object = "SnpAnnotationDataFrame"),
          function(object, ...) {
            getVariable(object, object@alleleBCol, ...)
          })

setMethod("getVariableNames",
          signature(object = "SnpAnnotationDataFrame"),
          function(object) {
            varLabels(object)
          })

setMethod("getAnnotation",
          signature(object = "SnpAnnotationDataFrame"),
          function(object) {
            pData(object)
          })

setMethod("getMetadata",
          signature(object = "SnpAnnotationDataFrame"),
          function(object) {
            varMetadata(object)
          })

setMethod("autosomeCode", "SnpAnnotationDataFrame",
          function(object) {
            object@autosomeCode
          })
   
setMethod("XchromCode", "SnpAnnotationDataFrame",
          function(object) {
            object@XchromCode
          })
          
setMethod("YchromCode", "SnpAnnotationDataFrame",
          function(object) {
            object@YchromCode
          })
              
setMethod("XYchromCode", "SnpAnnotationDataFrame",
          function(object) {
            object@XYchromCode
          })
              
setMethod("MchromCode", "SnpAnnotationDataFrame",
          function(object) {
            object@MchromCode
          })
          
