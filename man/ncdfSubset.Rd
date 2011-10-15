\name{ncdfSubset}
\alias{ncdfSubset}
\alias{ncdfSubsetCheck}
\title{Write a subset of data in a netCDF file to a new netCDF file}

\description{
  \code{ncdfSubset} takes a subset of data (snps and samples)
  from a netCDF file and write it to a new netCDF file.
  \code{ncdfSubsetCheck} checks that a netCDF file is the desired subset
  of another netCDF file.
}

\usage{
ncdfSubset(parent.ncdf, sub.ncdf,
           sample.include=NULL, snp.include=NULL,
           verbose=TRUE)

ncdfSubsetCheck(parent.ncdf, sub.ncdf,
                sample.include=NULL, snp.include=NULL,
                verbose=TRUE)
}

\arguments{
  \item{parent.ncdf}{Name of the parent netCDF file}
  \item{sub.ncdf}{Name of the subset netCDF file}
  \item{sample.include}{Vector of sampleIDs to include in \code{sub.ncdf}}
  \item{snp.include}{Vector of snpIDs to include in \code{sub.ncdf}}
  \item{verbose}{Logical value specifying whether to show progress information.}
}

\details{
  \code{ncdfSubset} can select a subset of snps for all samples by
  setting \code{snp.include}, a subset of samples for all snps by setting
  \code{sample.include}, or a subset of snps and samples with both arguments.
}

\author{Cathy Laurie, Stephanie Gogarten}

\seealso{\code{\link{ncdf}}, \code{\link{ncdfCreate}}, \code{\link{ncdfAddData}}
}

\examples{
ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(ncfile)
sample.sel <- getScanID(nc, index=1:10)
snp.sel <- getSnpID(nc, index=1:100)
close(nc)

subnc <- tempfile()
ncdfSubset(ncfile, subnc, sample.include=sample.sel, snp.include=snp.sel)
ncdfSubsetCheck(ncfile, subnc, sample.include=sample.sel, snp.include=snp.sel)
file.remove(subnc)
}

\keyword{manip}
