\name{NcdfReader}
\docType{class}

\alias{NcdfReader-class}
\alias{NcdfReader}
\alias{open,NcdfReader-method}
\alias{close,NcdfReader-method}
\alias{show,NcdfReader-method}
\alias{getDimensionNames}
\alias{getDimensionNames,NcdfReader-method}
\alias{getVariableNames}
\alias{getVariableNames,NcdfReader-method}
\alias{hasCoordVariable}
\alias{hasCoordVariable,NcdfReader-method}
\alias{hasVariable,NcdfReader-method}
\alias{getVariable,NcdfReader-method}
\alias{getAttribute}
\alias{getAttribute,NcdfReader-method}

\title{Class NcdfReader}

\description{
  The NcdfReader class is a wrapper for the \pkg{\link{ncdf}} library
  that provides an interface for reading NetCDF files.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{NcdfReader(filename)}:

      \code{filename} must be the path to a NetCDF file.

      The \code{NcdfReader} constructor creates and returns a
      NcdfReader instance pointing to this file.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a NcdfReader object.

  \describe{
    \item{}{
      \code{getVariable(object, varname, start, count)}: Returns the
      contents of the variable \code{varname}.     
      \itemize{
	\item \code{start} is a vector
	of integers indicating where to start reading values.  The length
	of this vector must equal the number of dimensions the variable has.
	If not specified, reading starts at the beginning of the file
	(1,1,...).
	\item \code{count} is a vector of integers indicating the count
	of values to read along each dimension.  The length of this
	vector must equal the number of dimensions the variable has. If not
	specified and the variable does NOT have an unlimited dimension, the
	entire variable is read. As a special case, the value "-1" indicates
	that all entries along that dimension should be read.
      }
      The result is a vector, matrix, or array, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.  If the variable is not found in the NetCDF
      file, returns \code{NULL}.
    }
    \item{}{
      \code{getVariableNames(object)}: Returns names of variables in the
      NetCDF file.
    }
    \item{}{
      \code{getDimensionNames(object, varname)}: Returns names of
      dimensions in the NetCDF file.  If \code{varname} is provided, returns
      dimension names for NetCDF variable \code{varname}.
    }
    \item{}{
      \code{getAttribute(object, attname, varname)}: Returns the
      attribute \code{attname} associated with the variable
      \code{varname}.  If \code{varname} is not specified, \code{attname}
      is assumed to be a global attribute.
      }
    \item{}{
      \code{hasCoordVariable(object, varname)}: Returns \code{TRUE} if
      \code{varname} is a coordinate variable (a variable with the same
      name as a dimension).
    }
    \item{}{
      \code{hasVariable(object, varname)}: Returns \code{TRUE} if
      \code{varname} is a variable in the NetCDF file (including
      coordinate variables).
    }
    \item{}{
      \code{open(object)}: Opens a connection to the NetCDF file.
    }
    \item{}{
      \code{close(object)}: Closes the NetCDF file connection.
    }
  }
}

\section{Standard Generic Methods}{
  In the code snippets below, \code{object} is a NcdfReader object.

  \describe{
    \item{}{
      \code{open(object)}: Opens a connection to a NetCDF file.
    }
    \item{}{
      \code{close(object)}: Closes the connection to a NetCDF file.
    }
  }
}

\author{Stephanie Gogarten}

\seealso{\pkg{\link{ncdf}}, \code{\link{NcdfGenotypeReader}},
  \code{\link{NcdfIntensityReader}}
}

\examples{
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfReader(file)

getDimensionNames(nc)
getVariableNames(nc)

hasVariable(nc, "genotype")
geno <- getVariable(nc, "genotype", start=c(1,1), count=c(10,10))

close(nc)
}

\keyword{methods}
\keyword{classes}
