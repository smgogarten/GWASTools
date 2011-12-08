\name{convertNcdfGds}
\alias{convertNcdfGds}
\alias{convertGdsNcdf}
\alias{checkNcdfGds}

\title{Convert between NetCDF and GDS format}

\description{
  \code{convertNcdfGds} converts a genotype NetCDF file to the GDS
  format used by the R packages \code{gdsfmt} and
  \code{SNPRelate}.

  \code{convertGdsNcdf} converts a GDS file to NetCDF format.

  \code{checkNcdfGds} checks whether a genotype NetCDF file and a GDS
  file contain identical data.
}

\usage{
convertNcdfGds(ncdf.filename, gds.filename,
  sample.annot = NULL, snp.annot = NULL, rsID.col = "rsID",
  alleleA.col = "alleleA", alleleB.col = "alleleB",
  zipflag = "zip.max", verbose = TRUE)

convertGdsNcdf(gds.filename, ncdf.filename, verbose = TRUE)

checkNcdfGds(ncdf.filename, gds.filename, verbose = TRUE)
}

\arguments{
  \item{ncdf.filename}{name of the NetCDF genotype file}
  \item{gds.filename}{name of the GDS file}
  \item{sample.annot}{a data.frame with sample annotation}
  \item{snp.annot}{a data.frame with SNP annotation}
  \item{rsID.col}{the name of the \code{snp.annot} column with rs ID}
  \item{alleleA.col}{the name of the \code{snp.annot} column with allele A}
  \item{alleleB.col}{the name of the \code{snp.annot} column with allele B}
  \item{zipflag}{the compression format for the GDS file, one of "",
    "ZIP", "ZIP.fast", "ZIP.default", or "ZIP.max"}
  \item{verbose}{whether to show progress information}
}

\details{
  These functions require that the package \code{gdsfmt} be
  installed.  \code{convertNcdfGds} is needed to convert the NetCDF
  genotype files used in this package to the format required by
  \code{SNPRelate} for Principal Component Analysis,
  Identity-by-Descent, and Linkage Disequlibrium calculations.
}

\value{
  \code{checkNcdfGds} returns \code{TRUE} if the NetCDF and GDS files
  contain identical data.  If the files differ, it will print a
  diagnostic message and return \code{FALSE}.
}

\author{Xiuwen Zheng}

\seealso{\code{gdsfmt}, \code{SNPRelate}, \code{\link{ncdf}},
  \code{\link{NcdfGenotypeReader}},
}

\examples{
library(GWASdata)
ncfile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")

data(illumina_snp_annot)
data(illumina_scan_annot)

gdsfile <- tempfile()
convertNcdfGds(ncfile, gdsfile, sample.annot=illumina_scan_annot,
  snp.annot=illumina_snp_annot, rsID.col="rsID",
  alleleA.col="allele.A", alleleB.col="allele.B")

checkNcdfGds(ncfile, gdsfile)

ncfile2 <- tempfile()
convertGdsNcdf(gdsfile, ncfile2)

file.remove(gdsfile, ncfile2)
}
\keyword{manip}