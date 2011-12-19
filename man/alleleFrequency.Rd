\name{alleleFrequency}
\alias{alleleFrequency}
\title{
Allelic frequency
}
\description{
Calculates the frequency of the A allele over the specifed scans.
}
\usage{
alleleFrequency(genoData, scan.exclude,
                verbose = TRUE)
}
\arguments{
  \item{genoData}{\code{\link{GenotypeData}} object.
}
  \item{scan.exclude}{Integer vector with IDs of scans to exclude.
}
  \item{verbose}{Logical value specifying whether to show progress information.
}
}
\details{
 Counts male heterozygotes on the X and Y chromosomes as missing values,
 and any genotype for females on the Y chromosome as missing values.
 A "sex" variable must be present in the scan
  annotation slot of \code{genoData}.
}
\value{
A matrix of allelic frequencies with snps as
rows and 3 columns ("M" for males, "F" for females, "all" for all scans).
}
\author{Cathy Laurie
}

\seealso{\code{\link{GenotypeData}}
}
\examples{
library(GWASdata)
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(file)

# need scan annotation with sex
data(affyScanADF)
genoData <-  GenotypeData(nc, scanAnnot=affyScanADF)

afreq <- alleleFrequency(genoData, scan.exclude=(affyScanADF$race != "CEU"))
close(genoData)
}
\keyword{manip}
