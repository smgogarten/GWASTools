\name{qualityScoreBySnp}
\alias{qualityScoreBySnp}

\title{ Mean and median quality score for SNPs }
\description{
This function calculates the mean and median quality score, over all scans with a non-missing genotype call, for each SNP.
}
\usage{
qualityScoreBySnp(intenData, genoData, scan.exclude = NULL,
                  block.size = 5000, verbose = TRUE)
}

\arguments{

  \item{intenData}{\code{\link{IntensityData}} object }
  \item{genoData}{\code{\link{GenotypeData}} object }
  \item{scan.exclude}{ An integer vector containing the id's of scans to be excluded. }
  \item{block.size}{Number of SNPs to be read from intenData and genoData at once.}
  \item{verbose}{Logical value specifying whether to show progress information.}

}



\details{
  \code{intenData} and \code{genoData} must have matching snpID
  and scanID.
}



\value{
  The function returns a matrix with the following columns:
  \item{mean.quality }{A vector of mean quality scores for each snp. }
  \item{median.quality }{A vector of median quality scores for each snp. }
}


\author{ Cathy Laurie }
\seealso{\code{\link{IntensityData}}, \code{\link{GenotypeData}}, \code{\link{qualityScoreByScan}} }


\examples{
qualfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
qual <- GdsIntensityReader(qualfile)
qualData <- IntensityData(qual)

genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
geno <- GdsGenotypeReader(genofile)
genoData <- GenotypeData(geno)

quality <- qualityScoreBySnp(qualData, genoData)
close(qualData)
close(genoData)
}


\keyword{univar }
