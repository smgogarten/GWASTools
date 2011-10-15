\name{getVariable}
\alias{getVariable}
\alias{getAnnotation}
\alias{getBAlleleFreq}
\alias{getChromosome}
\alias{getGenotype}
\alias{getLogRRatio}
\alias{getMetadata}
\alias{getPosition}
\alias{getQuality}
\alias{getScanID}
\alias{getScanVariable}
\alias{getScanVariableNames}
\alias{getSex}
\alias{getSnpID}
\alias{getSnpVariable}
\alias{getSnpVariableNames}
\alias{getX}
\alias{getY}
\alias{hasVariable}
\alias{hasBAlleleFreq}
\alias{hasLogRRatio}
\alias{hasQuality}
\alias{hasScanAnnotation}
\alias{hasScanVariable}
\alias{hasSex}
\alias{hasSnpAnnotation}
\alias{hasSnpVariable}
\alias{hasX}
\alias{hasY}
\alias{nsnp}
\alias{nscan}
\alias{writeAnnotation}
\alias{writeMetadata}
\alias{getQuery}
\alias{XchromCode}
\alias{YchromCode}
\alias{XYchromCode}
\alias{MchromCode}

\title{Accessors for variables in GenotypeData and IntensityData classes and
  their component classes}

\description{
  These generic functions provide access to variables associated with
  GWAS data cleaning.
}

\usage{
  getScanVariable(object, varname, ...)
  getScanID(object, ...)
  getSex(object, ...)
  getSnpVariable(object, varname, ...)
  getSnpID(object, ...)
  getChromosome(object, ...)
  getPosition(object, ...)
  getVariable(object, varname, ...)
  getGenotype(object, ...)
  getQuality(object, ...)
  getX(object, ...)
  getY(object, ...)
  getBAlleleFreq(object, ...)
  getLogRRatio(object, ...)
  getAnnotation(object, ...)
  getMetadata(object, ...)
  getQuery(object, statement, ...)

  hasScanAnnotation(object)
  hasScanVariable(object, varname)
  hasSex(object)
  hasSnpAnnotation(object)
  hasSnpVariable(object, varname)
  hasVariable(object, varname)
  hasQuality(object)
  hasX(object)
  hasY(object)
  hasBAlleleFreq(object)
  hasLogRRatio(object)

  nsnp(object)
  nscan(object)

  XchromCode(object)
  XYchromCode(object)
  YchromCode(object)
  MchromCode(object)

  writeAnnotation(object, value, ...)
  writeMetadata(object, value, ...)
}

\arguments{
  \item{object}{Object, possibly derived from or containing
    \code{\link{NcdfReader-class}},
    \code{\link{ScanAnnotationDataFrame-class}},
    \code{\link{SnpAnnotationDataFrame-class}},
    \code{\link{ScanAnnotationSQLite-class}}, or
    \code{\link{SnpAnnotationSQLite-class}}.}

  \item{varname}{Name of the variable (single character string, or a
    character vector for multiple variables).}

  \item{statement}{SQL statement to query
    \code{\link{ScanAnnotationSQLite-class}} or
    \code{\link{SnpAnnotationSQLite-class}} objects.}

  \item{value}{data.frame with annotation or metadata to write to
  \code{\link{ScanAnnotationSQLite-class}} or
  \code{\link{SnpAnnotationSQLite-class}} objects.}

  \item{...}{Additional arguments.}
}
  
\value{
  \code{get} methods return vectors or matrices of the requested
  variables (with the exception of \code{getQuery}, which returns a data
  frame).

  \code{has} methods return TRUE if the requested variable
  is present in \code{object}.

  \code{nsnp} and \code{nscan} return the number of SNPs and scans in
  the object, repectively.

  \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
  \code{MchromCode} return the integer chromosome codes associated with X,
  pseudoautosomal, Y, and mitochondrial SNPs.
}

\author{Stephanie Gogarten}

\seealso{\code{\link{ScanAnnotationDataFrame-class}},
  \code{\link{SnpAnnotationDataFrame-class}},
  \code{\link{ScanAnnotationSQLite-class}},
  \code{\link{SnpAnnotationSQLite-class}},
  \code{\link{NcdfReader-class}},
  \code{\link{NcdfGenotypeReader-class}},
  \code{\link{NcdfIntensityReader-class}},
  \code{\link{GenotypeData-class}}, \code{\link{IntensityData-class}}
}

\keyword{methods}
