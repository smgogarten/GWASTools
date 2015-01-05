# Methods for GdsReader

# constructor
GdsReader <- function(filename) {
  if (missing(filename)) stop("filename is required")
  if (is(filename, 'gds.class')) {
    handler <- filename
    filename <- handler$filename
  } else {
    if (!file.exists(filename)) stop("Error in opening file ", filename, ": no such file or directory")
    handler <- openfn.gds(filename)
  }
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
            print(object@handler)
          })

setMethod("getNodeDescription",
          signature(object="GdsReader"),
          function(object, varname) {
            objdesp.gdsn(index.gdsn(object@handler, varname))
          })

setMethod("getDimension",
          signature(object="GdsReader"),
          function(object, varname) {
            getNodeDescription(object, varname)$dim
          })

setMethod("getVariableNames",
          signature(object="GdsReader"),
          function(object) {
            vars <- ls.gdsn(object@handler)
            # number of child nodes
            n.child <- sapply(vars, function(x) cnt.gdsn(index.gdsn(object@handler, x)))
            folders <- vars[n.child > 0]
            if (length(folders) > 0) {
                varf <- unlist(lapply(folders, function(x)
                    paste(x, ls.gdsn(index.gdsn(object@handler, x)), sep="/")))
                vars <- c(setdiff(vars, folders), varf)
            }
            vars
          })

setMethod("hasVariable",
          signature(object="GdsReader"),
          function(object, varname) {
            varname %in% getVariableNames(object)
          })

setMethod("getVariable",
          signature(object="GdsReader"),
          function(object, varname, sel=NULL, ...) {

            # check that variable exists
            if (!hasVariable(object, varname)) {
              warning(paste(varname, "not found"))
              return(NULL)
            }

            # get variable from gds
            node <- index.gdsn(object@handler, varname)
            if (is.null(sel)) {
                var <- read.gdsn(node, ...)
            } else {
                var <- readex.gdsn(node, sel, ...)
            }

            # set missing value to NA
            missVal <- getAttribute(object, "missing.value", varname)
            if (!is.null(missVal)) {
              var[var == missVal] <- NA
            }

            return(var)
          })

setMethod("getAttribute",
          signature(object="GdsReader"),
          function(object, attname, varname) {
            get.attr.gdsn(index.gdsn(object@handler, varname))[[attname]]
          })
