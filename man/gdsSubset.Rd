\name{gdsSubset}
\alias{gdsSubset}
\alias{gdsSubsetCheck}
\title{Write a subset of data in a GDS file to a new GDS file}

\description{
  \code{gdsSubset} takes a subset of data (snps and samples)
  from a GDS file and write it to a new GDS file.
  \code{gdsSubsetCheck} checks that a GDS file is the desired subset
  of another GDS file.
}

\usage{
gdsSubset(parent.gds, sub.gds,
          sample.include=NULL, snp.include=NULL,
          sub.storage=NULL,
          zipflag="ZIP.max",
          verbose=TRUE)

gdsSubsetCheck(parent.gds, sub.gds,
               sample.include=NULL, snp.include=NULL,
               sub.storage=NULL,
               verbose=TRUE)
}

\arguments{
  \item{parent.gds}{Name of the parent GDS file}
  \item{sub.gds}{Name of the subset GDS file}
  \item{sample.include}{Vector of sampleIDs to include in \code{sub.gds}}
  \item{snp.include}{Vector of snpIDs to include in \code{sub.gds}}
  \item{sub.storage}{storage type for the subset file; defaults to original storage type}
  \item{zipflag}{compression for GDS file, can be "", "ZIP", "ZIP.fast", "ZIP.default", "ZIP.max"}
  \item{verbose}{Logical value specifying whether to show progress information.}

}

\details{
  \code{gdsSubset} can select a subset of snps for all samples by
  setting \code{snp.include}, a subset of samples for all snps by setting
  \code{sample.include}, or a subset of snps and samples with both arguments.
}

\author{Adrienne Stilp}

\seealso{\pkg{\link{gdsfmt}}, \code{convertNcdfGds}}

\examples{
ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
gdsfile <- tempfile()
convertNcdfGds(ncdf.filename=ncfile, gds.filename=gdsfile)

gds <- GdsGenotypeReader(gdsfile)
sample.sel <- getScanID(gds, index=1:10)
snp.sel <- getSnpID(gds, index=1:100)
close(gds)

subfile <- tempfile()
gdsSubset(gdsfile, subfile, sample.include=sample.sel, snp.include=snp.sel)
gdsSubsetCheck(gdsfile, subfile, sample.include=sample.sel, snp.include=snp.sel)

file.remove(c(gdsfile, subfile))
}

\keyword{manip}
