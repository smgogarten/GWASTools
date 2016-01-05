# Methods for NcdfReader

# constructor
NcdfReader <- function(filename) {
  if (missing(filename)) stop("filename is required")
  if (!file.exists(filename)) stop("Error in opening file ", filename, ": no such file or directory")
  handler <- nc_open(filename, readunlim=FALSE)
  new("NcdfReader", filename=filename, handler=handler)
}

setValidity("NcdfReader",
            function(object) {
              if (!is.character(object@filename) ||
                  length(object@filename) != 1 ||
                  is.na(object@filename))
                return("'filename' slot must be a single string")
              TRUE
            })

setMethod("open",
    signature(con = "NcdfReader"),
    function (con, ...) {
      con@handler <- nc_open(con@filename, readunlim=FALSE, ...)
    })

setMethod("close",
    signature(con = "NcdfReader"),
    function (con, ...) {
      x <- nc_close(con@handler, ...)
    })

setMethod("show",
          signature(object="NcdfReader"),
          function(object) {
            print(object@handler)
          })

setMethod("getDimension",
          signature(object="NcdfReader"),
          function(object, varname) {
            sapply(object@handler$var[[varname]]$dim, function(x) x$len)
          })

# get dimension names
# if varname is missing, returns dimension names for netcdf object
setMethod("getDimensionNames",
          signature(object="NcdfReader"),
          function(object, varname) {
            if (missing(varname)) {
              names(object@handler$dim)
            } else {
              sapply(object@handler$var[[varname]]$dim, function(x) x$name)
            }
          })

setMethod("getVariableNames",
          signature(object="NcdfReader"),
          function(object) {
            names(object@handler$var)
          })

# returns TRUE if varname is a coordinate variable (variable with same
# name as a dimension)
setMethod("hasCoordVariable",
          signature(object="NcdfReader"),
          function(object, varname) {
            isDimension <- varname %in% getDimensionNames(object)
            if (isDimension) {
              object@handler$dim[[varname]]$create_dimvar
            } else {
              isDimension
            }
          })

# include both regular and coordinate variables
setMethod("hasVariable",
          signature(object="NcdfReader"),
          function(object, varname) {
            varname %in% getVariableNames(object) |
            hasCoordVariable(object, varname)
          })


setMethod("getVariable",
          signature(object="NcdfReader"),
          function(object, varname, start, count, drop=TRUE, ...) {

            # check that variable exists
            if (!hasVariable(object, varname)) {
              warning(paste(varname, "not found"))
              return(NULL)
            }

            # if start and count not specified, return all elements
            if (missing(start)) start <- NA
            if (missing(count)) count <- NA

            # get variable from netcdf
            var <- ncvar_get(object@handler, varname, start, count)

            # 1D variables are returned as arrays - convert to vector
            if (is(var, "array") & length(dim(var)) == 1) {
              if (drop) {
                var <- as.vector(var)
              } else {
                if (!all(is.na(count)) & length(count) == 2) {
                  if (count[1] == 1) {
                    var <- matrix(var, nrow=1)
                  } else if (count[2] == 1) {
                    var <- matrix(var, ncol=1)
                  }
                }
              }
            }

            return(var)
          })

setMethod("getAttribute",
          signature(object="NcdfReader"),
          function(object, attname, varname) {
            if (missing(varname)) {
              varname <- 0
            }
            ncatt_get(object@handler, varname, attname)$value
          })
