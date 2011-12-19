\name{qualityScoreByScan}
\alias{qualityScoreByScan}

\title{ Mean and median quality score for scans }
\description{
This function calculates the mean and median quality score, over all SNPs with a non-missing genotype call, for each scan.
}
\usage{
qualityScoreByScan(intenData, genoData,
                   snp.exclude = NULL,
                   verbose = TRUE)
}

\arguments{

  \item{intenData}{\code{\link{IntensityData}} object }
  \item{genoData}{\code{\link{GenotypeData}} object }
  \item{snp.exclude}{ An integer vector containing the id's of SNPs to be excluded. }
  \item{verbose}{Logical value specifying whether to show progress information.}

}



\details{
  \code{intenData} and \code{genoData} must have matching snpID
  and scanID.
  Y chromosome SNPs are excluded for females. 
  A "sex" variable must be present in the scan
  annotation slot of \code{intenData} or \code{genoData}.
}



\value{
  The function returns a matrix with the following columns:
  \item{mean.quality}{A vector of mean quality scores for each scan }
  \item{median.quality}{A vector of median quality scores for each scan. }
}


\author{ Cathy Laurie }

\seealso{\code{\link{IntensityData}}, \code{\link{GenotypeData}}, \code{\link{qualityScoreBySnp}}  }


\examples{
library(GWASdata)
qualfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
qualnc <- NcdfIntensityReader(qualfile)
# need scan annotation with sex
data(affyScanADF)
qualData <- IntensityData(qualnc, scanAnnot=affyScanADF)

genofile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
genonc <- NcdfGenotypeReader(genofile)
genoData <- GenotypeData(genonc, scanAnnot=affyScanADF)

quality <- qualityScoreByScan(qualData, genoData)
close(qualData)
close(genoData)
}


\keyword{univar }
