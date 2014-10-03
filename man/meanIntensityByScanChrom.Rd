\name{meanIntensityByScanChrom}
\alias{meanIntensityByScanChrom}
\title{Calculate Means \& Standard Deviations of Intensities}
\description{
Function to calculate the mean and standard deviation of the intensity for each chromosome for each scan.
}


\usage{
meanIntensityByScanChrom(intenData, vars = c("X", "Y"), 
                         snp.exclude = NULL, verbose = TRUE)
}


\arguments{
  \item{intenData}{\code{\link{IntensityData}} object.  Chromosomes are
  expected to be in contiguous blocks.}
  \item{vars}{Character vector with the names of one or two intensity
  variables.}
  \item{snp.exclude}{ An integer vector containing SNPs to be excluded. }
  \item{verbose}{Logical value specifying whether to show progress information.}
}


\details{
  The names of two intensity variables in \code{intenData} may be supplied.
  If two variables are given, the mean of their sum is computed as
  well.  The default is to compute the mean and standard deviation for X
  and Y intensity.
}


\value{
A list with two components for each variable in "vars": 'mean.var' and
'sd.var'.  If two variables are given, the first two elements of the
list will be mean and sd for the sum of the intensity variables:  
  \item{mean.intensity }{A matrix with one row per scan and one column
    per chromosome containing the means of the summed intensity values for each scan and chromosome. }
  \item{sd.intensity }{A matrix with one row per scan and one column
    per chromosome containing the standard deviations of the summed intensity values for each scan and chromosome. }
  \item{mean.var }{A matrix with one row per scan and one column per chromosome containing the means of the intensity values for each scan and chromosome. }
  \item{sd.var }{A matrix with one row per scan and one column per chromosome containing the standard deviations of the intensity values for each scan and chromosome. }
}
\author{ Cathy Laurie }

\seealso{\code{\link{IntensityData}}, \code{\link{mean}}, \code{\link{sd}}}
\examples{
file <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
gds <- GdsIntensityReader(file)
intenData <- IntensityData(gds)

meanInten <- meanIntensityByScanChrom(intenData)
close(intenData)
}
\keyword{univar}
