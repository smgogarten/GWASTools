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

Samples with missing sex are included in the allele counts for "all" and "MAF" for autosomes, but not for sex chromosomes.
}
\value{
A matrix with a row for each SNP.
Columns "M" for males, "F" for females, and "all" for all scans give frequencies of the A allele.
Sample size for males, females, and all is returned as "n.M", "n.F", and "n", respectively.
"MAF" is the minor allele frequency over all scans.
}
\author{Cathy Laurie, Stephanie Gogarten
}

\seealso{\code{\link{GenotypeData}}
}
\examples{
library(GWASdata)
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)

# need scan annotation with sex
data(illuminaScanADF)
genoData <-  GenotypeData(gds, scanAnnot=illuminaScanADF)

afreq <- alleleFrequency(genoData, scan.exclude=(illuminaScanADF$race != "CEU"))
close(genoData)
}
\keyword{manip}
