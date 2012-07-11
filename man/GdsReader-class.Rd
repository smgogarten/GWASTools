\name{GdsReader}
\docType{class}

\alias{GdsReader-class}
\alias{GdsReader}
\alias{open,GdsReader-method}
\alias{close,GdsReader-method}
\alias{show,GdsReader-method}
\alias{getDimension}
\alias{getDimension,GdsReader-method}
\alias{getVariableNames,GdsReader-method}
\alias{hasVariable,GdsReader-method}
\alias{getVariable,GdsReader-method}
\alias{getAttribute,GdsReader-method}

\title{Class GdsReader}

\description{
  The GdsReader class provides an interface for reading GDS files.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GdsReader(filename)}:

      \code{filename} must be the path to a GDS file.

      The \code{GdsReader} constructor creates and returns a
      GdsReader instance pointing to this file.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a GdsReader object.

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
      of dimensions in the returned values.  If the variable is not found in the GDS
      file, returns \code{NULL}.
    }
    \item{}{
      \code{getVariableNames(object)}: Returns names of variables in the
      GDS file.
    }
    \item{}{
      \code{getDimension(object, varname)}: Returns dimension
      for GDS variable \code{varname}.
    }
    \item{}{
      \code{getAttribute(object, attname, varname)}: Returns the
      attribute \code{attname} associated with the variable
      \code{varname}.
    }
    \item{}{
      \code{hasVariable(object, varname)}: Returns \code{TRUE} if
      \code{varname} is a variable in the GDS file.
    }
  }
}

\section{Standard Generic Methods}{
  In the code snippets below, \code{object} is a GdsReader object.

  \describe{
    \item{}{
      \code{open(object)}: Opens a connection to a GDS file.
    }
    \item{}{
      \code{close(object)}: Closes the connection to a GDS file.
    }
    \item{}{
      \code{show(object)}: Summarizes the contents of a GDS file.
    }
  }
}

\author{Stephanie Gogarten}

\seealso{\pkg{\link{gdsfmt}}
  %\code{\link{GdsGenotypeReader}}
}

\examples{
library(SNPRelate)
gds <- GdsReader(snpgdsExampleFileName())

getVariableNames(gds)

hasVariable(gds, "genotype")
geno <- getVariable(gds, "genotype", start=c(1,1), count=c(10,10))

close(gds)
}

\keyword{methods}
\keyword{classes}
