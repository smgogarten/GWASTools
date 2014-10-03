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
          block.size=5000,
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
  \item{block.size}{for GDS files stored with scan,snp dimensions, the number of SNPs to read from the parent file at a time. Ignored for snp,scan dimensions.}
  \item{verbose}{Logical value specifying whether to show progress information.}

}

\details{
  \code{gdsSubset} can select a subset of snps for all samples by setting \code{snp.include}, a subset of samples for all snps by setting \code{sample.include}, or a subset of snps and samples with both arguments.
  The GDS nodes "snp.id", "snp.position", "snp.chromosome", and "sample.id" are copied, as well as any 2-dimensional nodes. Other nodes are not copied.
  The attributes of the 2-dimensional nodes are also copied to the subset file.
  If sub.storage is specified, the subset gds file will have a different storage mode for any 2-dimensional array.
  In the special case where the 2-dimensional node has an attribute named \code{"missing.value"} and the sub.storage type is \code{"bit2"}, the missing.value attribute for the subset node is automatically set to 3.
  At this point, no checking is done to ensure that the values will be properly stored with a different storage type, but \code{\link{gdsSubsetCheck}} will return an error if the values do not match.
  If the nodes in the GDS file are stored with scan,snp dimensions, then \code{block.size} allows you to loop over a block of SNPs at a time.
  If the nodes are stored with snp,scan dimensions, then the function simply loops over samples, one at a time.
  
  \code{gdsSubsetCheck} checks that a subset GDS file has the expected SNPs and samples of the parent file. It also checks that attributes were similarly copied, except for the above-mentioned special case of \code{missing.value} for \code{sub.storage="bit2"}.
}

\author{Adrienne Stilp}

\seealso{\pkg{\link{gdsfmt}}, \code{\link{createDataFile}}}

\examples{
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(gdsfile)
sample.sel <- getScanID(gds, index=1:10)
snp.sel <- getSnpID(gds, index=1:100)
close(gds)

subfile <- tempfile()
gdsSubset(gdsfile, subfile, sample.include=sample.sel, snp.include=snp.sel)
gdsSubsetCheck(gdsfile, subfile, sample.include=sample.sel, snp.include=snp.sel)

file.remove(subfile)
}

\keyword{manip}
