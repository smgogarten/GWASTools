\name{asSnpMatrix}
\alias{asSnpMatrix}
\title{Utilities for snpStats}
\description{\code{asSnpMatrix} converts a \code{\link{GenotypeData}}
  object to a \code{\link[snpStats]{SnpMatrix-class}} object.}
\usage{
asSnpMatrix(genoData, snpNames="snpID", scanNames="scanID",
            snp=c(1,-1), scan=c(1,-1))
}
\arguments{
  \item{genoData}{A \code{\link{GenotypeData}} object.}
  \item{snpNames}{The name of the SNP variable in \code{genoData} to use
    as the column (SNP) names in the \code{\link[snpStats]{SnpMatrix-class}} object.}
  \item{scanNames}{The name of the scan variable in \code{genoData} to use
    as the row (scan) names in the \code{\link[snpStats]{SnpMatrix-class}} object.}
  \item{snp}{An integer vector of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that all SNPs should be read.}
  \item{scan}{An integer vector of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that all scans should be read.}
}
\details{
  The default is to extract all SNPs and scans from \code{genoData}, but
  for a large dataset this may exceed R's memory limit.  Alternatively,
  \code{snp} and \code{scan} may be used to specify \code{(start,
    count)} of SNPs and scans to extract from \code{genoData}.
  
  In the \code{SnpMatrix} object, genotypes are stored as 0 = missing, 1
  = "A/A", 2=  "A/B" or "B/A", and 3 = "B/B".  (In a \code{GenotypeData}
  object, 0 = "B/B", 1 = "A/B" or "B/A", and 2 = "A/A".)
  Columns are SNPs with
  names \code{snpNames} and rows are scans with names \code{scanNames}
  (the transpose of the \code{GenotypeData} object).
}
\value{
  A \code{\link[snpStats]{SnpMatrix-class}} object.
}
\author{Stephanie Gogarten}
\seealso{\code{\link[snpStats]{SnpMatrix-class}}, \code{\link{GenotypeData}}}
\examples{
library(snpStats)
library(GWASdata)
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(file)
data(affySnpADF)
data(affyScanADF)
genoData <-  GenotypeData(nc, snpAnnot=affySnpADF, scanAnnot=affyScanADF)
snpmat <- asSnpMatrix(genoData, snpNames="rsID", scanNames="scanID")
snpmat
as(snpmat[1:5, 1:5], "character")
summary(snpmat)

# only chromosome 21
chr <- getChromosome(genoData)
c21 <- which(chr == 21)
snpmat <- asSnpMatrix(genoData, snpNames="rsID", scanNames="scanID",
                      snp=c(c21[1], length(c21)))
snpmat
}
\keyword{manip}
