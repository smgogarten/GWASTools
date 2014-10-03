\name{missingGenotypeBySnpSex}
\alias{missingGenotypeBySnpSex}
\title{Missing Counts per SNP by Sex}
\description{
For all SNPs for each sex tabulates missing SNP counts, allele counts and heterozygous counts.  
}
\usage{
missingGenotypeBySnpSex(genoData, scan.exclude = NULL, 
                        verbose = TRUE)
}


\arguments{
  \item{genoData}{\code{\link{GenotypeData}} object.}
  \item{scan.exclude}{A vector containing the scan numbers of scans that are to be excluded from the total scan list. }
  \item{verbose}{Logical value specifying whether to show progress information.}
}


\details{
  This function calculates the fraction of missing genotypes for males
  and females for each SNP given in \code{genoData}.
  A "sex" variable must be present in the scan
  annotation slot of \code{genoData}.
}


\value{
This function returns a list with three components: "missing.counts,"
"scans.per.sex," and "missing.fraction."
  \item{missing.counts}{A matrix with one row per SNP and one column per sex containing the number of missing SNP counts for males and females, respectively. }
  \item{scans.per.sex}{A vector containing the number of males and
    females respectively.}
  \item{missing.fraction}{A vector containing the fraction of missing
    counts for each SNP, with females excluded for the Y chromosome.}
}
\author{ Cathy Laurie, Stephanie Gogarten }
\seealso{\code{\link{GenotypeData}}, \code{\link{missingGenotypeByScanChrom}}}
\examples{
library(GWASdata)
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)

# need scan annotation with sex
data(illuminaScanADF)
genoData <-  GenotypeData(gds, scanAnnot=illuminaScanADF)

missingRate <- missingGenotypeBySnpSex(genoData)
close(genoData)
}
\keyword{univar}