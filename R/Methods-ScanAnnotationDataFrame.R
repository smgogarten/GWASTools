# Methods for ScanAnnotationDataFrame

# constructor
ScanAnnotationDataFrame <- function(data, metadata, ...) {
  if (missing(metadata)) {
    new("ScanAnnotationDataFrame", data=data, dimLabels=c("scans", "variables"), ...)
  } else {
    new("ScanAnnotationDataFrame", data=data, varMetadata=metadata,
        dimLabels=c("scans", "variables"), ...)
  }
}

# validity
setValidity("ScanAnnotationDataFrame",
          function(object) {
            # check for required columns
            if (!hasVariable(object, object@idCol)) {
              return(paste("missing required column", object@idCol))
            }

            # check format of required columns
            # scanID should be a unique vector
            scanID <- getScanID(object)
            if (length(scanID) != length(unique(scanID))) {
              return(paste(object@idCol, "must be a unique vector"))
            }
            # sex should be M/F
            if (hasSex(object)) {
              sex <- getSex(object)
              if (!all(sex %in% c("M","F") | is.na(sex))) {
                return(paste(object@sexCol, "should have values M/F"))
              }
            }
            TRUE
          })
     
setMethod("hasVariable",
          signature(object = "ScanAnnotationDataFrame"),
          function(object, varname) {
            varname %in% varLabels(object)
          })
            
setMethod("getVariable",
          signature(object = "ScanAnnotationDataFrame"),
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

setMethod("getScanID",
          signature(object = "ScanAnnotationDataFrame"),
          function(object, ...) {
            getVariable(object, object@idCol, ...)
          })
  
setMethod("hasSex",
          signature(object = "ScanAnnotationDataFrame"),
          function(object) {
            hasVariable(object, object@sexCol)
          })
        
setMethod("getSex",
          signature(object = "ScanAnnotationDataFrame"),
          function(object, ...) {
            getVariable(object, object@sexCol, ...)
          })

setMethod("getVariableNames",
          signature(object = "ScanAnnotationDataFrame"),
          function(object) {
            varLabels(object)
          })

setMethod("getAnnotation",
          signature(object = "ScanAnnotationDataFrame"),
          function(object) {
            pData(object)
          })

setMethod("getMetadata",
          signature(object = "ScanAnnotationDataFrame"),
          function(object) {
            varMetadata(object)
          })
