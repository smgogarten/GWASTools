\name{ncdfCreate}
\alias{ncdfCreate}


\title{
Write genotypic calls and/or associated metrics to a netCDF file - deprecated
}

\description{
Deprecated - use \code{\link{createDataFile}}

The function creates a shell netCDF file to which data can subsequently written. 
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usage{
ncdfCreate(snp.annotation, ncdf.filename, variables = "genotype", 
           n.samples = 10, precision = "double",
           array.name = NULL, genome.build = NULL)
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\arguments{

\item{snp.annotation}{Snp annotation dataframe with columns "snpID",
  "chromosome", and "position".  snpID should be a unique integer
  vector, sorted with
  respect to chromosome and position.}

\item{ncdf.filename}{The name of the genotype netCDF file to create}
 
\item{variables}{A character vector containing the names of the
  variables to create (must be one or more of \code{c("genotype", "quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio")})}

\item{n.samples}{The number of samples that will be in the netcdf file.}

\item{precision}{A character value indicating whether floating point
  numbers should be stored as "double" or "single" precision.}

\item{array.name}{Name of the array, to be stored as an attribute in the
  netCDF file.}

\item{genome.build}{Genome build used in determining chromosome and
  position, to be stored as an attribute in the netCDF file.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\details{
The function creates a shell netCDF file to which data can subsequently written. 
  
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\author{
Cathy Laurie
}
\seealso{\code{\link{ncdf}}, \code{\link{ncdfAddData}}, \code{\link{ncdfSubset}}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
\dontrun{
library(GWASdata)
data(affy_snp_annot)
ncfile <- tempfile()
ncdfCreate(affy_snp_annot, ncfile, variables="genotype", n.samples=5)
file.remove(ncfile)
}
}

\keyword{manip}
