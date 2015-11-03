\name{convertNcdfGds}
\alias{convertNcdfGds}
\alias{convertGdsNcdf}
\alias{checkNcdfGds}

\title{Convert between NetCDF and GDS format}

\description{
  \code{convertNcdfGds} converts a NetCDF file to GDS
  format.

  \code{convertGdsNcdf} converts a GDS file to NetCDF format.

  \code{checkNcdfGds} checks whether a genotype NetCDF file and a GDS
  file contain identical data.
}

\usage{
convertNcdfGds(ncdf.filename, gds.filename, snp.annot = NULL,
  precision = "single", compress = "ZIP_RA", verbose = TRUE)

convertGdsNcdf(gds.filename, ncdf.filename,
  precision = "single", verbose = TRUE)

checkNcdfGds(ncdf.filename, gds.filename, verbose = TRUE)
}

\arguments{
  \item{ncdf.filename}{name of the NetCDF file}
  \item{gds.filename}{name of the GDS file}
  \item{snp.annot}{a \code{\link{SnpAnnotationDataFrame}} with SNP
    annotation.  The column named "snpName" will be written to
  "snp.rs.id" in the GDS file.}

\item{precision}{A character value indicating whether floating point
  numbers should be stored as "double" or "single" precision.}

\item{compress}{The compression level for variables in a GDS file (see
\code{\link{add.gdsn}} for options.}
  \item{verbose}{whether to show progress information}
}

\details{
  \code{convertNcdfGds} assumes any variables other than "sampleID",
  "chromosome", and "position" have dimensions SNP x sample.

  If \code{snp.annot} has columns "rsID", "alleleA", "alleleB", these
  will be stored in the GDS file as "snp.rs.id" and "snp.allele" (the
  latter has the format "A/B").

  Chromosome codes from \code{snp.annot} (for autosomes, X, Y, etc.)
  will be stored in the GDS file.

  \code{convertGdsNcdf} assumes any variables not starting with "snp" or
  "sample" have dimensions SNP x sample.
}

\value{
  \code{checkNcdfGds} returns \code{TRUE} if the NetCDF and GDS files
  contain identical data.  If the files differ, it will print a
  diagnostic message and return \code{FALSE}.
}

\author{Xiuwen Zheng, Stephanie Gogarten}

\seealso{\pkg{\link[gdsfmt]{gdsfmt}},  \pkg{\link{ncdf}}
}

\examples{
library(GWASdata)
ncfile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")

data(illuminaSnpADF)

gdsfile <- tempfile()
convertNcdfGds(ncfile, gdsfile, snp.annot=illuminaSnpADF)

checkNcdfGds(ncfile, gdsfile)

ncfile2 <- tempfile()
convertGdsNcdf(gdsfile, ncfile2)

file.remove(gdsfile, ncfile2)
}
\keyword{manip}
