\name{ncdfSetMissingGenotypes}
\alias{ncdfSetMissingGenotypes}
\title{Write a new netCDF file, setting certain SNPs to missing}

\description{
  \code{ncdfSetMissingGenotypes} copies an existing netCDF genotype file to a new
  one, setting SNPs in specified regions to missing.
}

\usage{
ncdfSetMissingGenotypes(parent.ncdf, new.ncdf, regions,
                        sample.include=NULL, verbose=TRUE)
}

\arguments{
  \item{parent.ncdf}{Name of the parent netCDF file}
  \item{new.ncdf}{Name of the new netCDF file}
  \item{regions}{Data.frame of chromosome regions with columns
  \code{"scanID", "chromosome", "left.base", "right.base",
  "whole.chrom"}.}
  \item{sample.include}{Vector of sampleIDs to include in \code{new.ncdf}}
  \item{verbose}{Logical value specifying whether to show progress information.}
}

\details{
  \code{ncdfSetMissingGenotypes} removes chromosome regions by setting
  SNPs that fall within the anomaly regions to \code{NA} (i.e., the missing value
  in the netCDF file).  Optionally, entire samples may be excluded from
  the netCDF file as well: if the \code{sample.include} argument is
  given, only the scanIDs in this vector will be written to the new
  file, so the sample dimension will be \code{length(sample.include)}.

  For regions with \code{whole.chrom=TRUE}, the entire chromosome will
  be set to \code{NA} for that sample.  For other regions, only the
  region between \code{left.base} and \code{right.base} will be set to \code{NA}.
}

\author{Stephanie Gogarten}

\seealso{\code{\link{ncdfSubset}}, \code{\link{anomSegStats}} for
  chromosome anomaly regions
}

\examples{
ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(ncfile)
sample.sel <- getScanID(nc, index=1:10)
close(nc)

regions <- data.frame("scanID"=sample.sel[1:3], "chromosome"=c(21,22,23),
  "left.base"=c(14000000, 30000000, NA), "right.base"=c(28000000, 450000000, NA),
  whole.chrom=c(FALSE, FALSE, TRUE))

newnc <- tempfile()
ncdfSetMissingGenotypes(ncfile, newnc, regions, sample.include=sample.sel)
file.remove(newnc)
}

\keyword{manip}
