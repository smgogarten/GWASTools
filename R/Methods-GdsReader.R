# Methods for GdsReader

# constructor
GdsReader <- function(filename) {
  if (missing(filename)) stop("filename is required")
  if (!file.exists(filename)) stop("Error in opening file ", filename, ": no such file or directory")
  handler <- openfn.gds(filename)
  new("GdsReader", filename=filename, handler=handler)
}

setValidity("GdsReader",
            function(object) {              
              if (!is.character(object@filename) ||
                  length(object@filename) != 1 ||
                  is.na(object@filename))
                return("'filename' slot must be a single string")
              TRUE
            })

setMethod("open",
    signature(con = "GdsReader"),
    function (con) {
      con@handler <- openfn.gds(con@filename)
    })

setMethod("close",
    signature(con = "GdsReader"),
    function (con) {
      x <- closefn.gds(con@handler)
    })

setMethod("show", 
          signature(object="GdsReader"),
          function(object) {
            show(object@handler)
          })

setMethod("getDimension",
          signature(object="GdsReader"),
          function(object, varname) {
            objdesp.gdsn(index.gdsn(object@handler, varname))$dim
          })

setMethod("getVariableNames",
          signature(object="GdsReader"),
          function(object) {
            ls.gdsn(object@handler)
          })

setMethod("hasVariable",
          signature(object="GdsReader"),
          function(object, varname) {
            varname %in% getVariableNames(object)
          })

setMethod("getVariable",
          signature(object="GdsReader"),
          function(object, varname, ...) {
            
            # check that variable exists
            if (!hasVariable(object, varname)) {
              warning(paste(varname, "not found"))
              return(NULL)
            }

            # get variable from gds
            var <- read.gdsn(index.gdsn(object@handler, varname), ...)

            return(var)
          })

setMethod("getAttribute",
          signature(object="GdsReader"),      
          function(object, attname, varname) {
            get.attr.gdsn(index.gdsn(object@handler, varname)[[attname]])
          })
