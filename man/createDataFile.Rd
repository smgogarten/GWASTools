\name{createDataFile}
\alias{createDataFile}


\title{
Write genotypic calls and/or associated metrics to a GDS or netCDF file.
}

\description{
The function creates a shell GDS or netCDF file to which data can subsequently written.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usage{
createDataFile(snp.annotation, filename, file.type=c("gds", "ncdf"),
               variables="genotype", n.samples=10, precision="single",
               compress="ZIP.max", array.name = NULL, genome.build = NULL)
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\arguments{

\item{snp.annotation}{Snp annotation dataframe with columns "snpID",
  "chromosome", and "position".  snpID should be a unique integer
  vector, sorted with
  respect to chromosome and position.  If \code{file.type="gds"},
  optional columns named "rsID", "alleleA", and "alleleB" will be
  written if present.}

\item{filename}{The name of the genotype GDS or netCDF file to create}

\item{file.type}{The type of file to create ("gds" or "ncdf")}

\item{variables}{A character vector containing the names of the
  variables to create (must be one or more of \code{c("genotype", "quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio")})}

\item{n.samples}{The number of samples that will be in the netCDF file.
Not needed for a GDS file, as the file is expanded as each sample is written.}

\item{precision}{A character value indicating whether floating point
  numbers should be stored as "double" or "single" precision.}

\item{compress}{The compression level for variables in a GDS file (see
\code{\link{add.gdsn}} for options.}

\item{array.name}{Name of the array, to be stored as an attribute in the
  netCDF file.}

\item{genome.build}{Genome build used in determining chromosome and
  position, to be stored as an attribute in the netCDF file.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\details{
The function creates a shell GDS or netCDF file to which data can subsequently written.

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\author{
Stephanie Gogarten, Cathy Laurie
}
\seealso{\code{\link{ncdf}}, \code{\link{ncdfAddData}}, \code{\link{ncdfSubset}}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
library(GWASdata)
data(affy_snp_annot)
gdsfile <- tempfile()
createDataFile(affy_snp_annot, gdsfile, file.type="gds", variables="genotype", n.samples=5)
file.remove(gdsfile)
}

\keyword{manip}
