\name{intensityOutliersPlot}
\alias{intensityOutliersPlot}
\title{Plot mean intensity and highlight outliers}
\description{
  \code{intensityOutliersPlot} is a function to plot mean
  intensity for chromosome i vs mean of intensities for autosomes
  (excluding i) and highlight outliers
}
\usage{
intensityOutliersPlot(mean.intensities, sex, outliers,	
	              sep = FALSE, label, ...)
}
\arguments{
  \item{mean.intensities}{scan x chromosome matrix of mean intensities}
  \item{sex}{vector with values of "M" or "F" corresponding to scans in the rows of \code{mean.intensities}}
  \item{outliers}{list of outliers, each member corresponds to a chromosome (member "X" is itself a list of female and male outliers)}
  \item{sep}{plot outliers within a chromosome separately (TRUE) or together (FALSE)}
  \item{label}{list of plot labels (to be positioned below X axis)
    corresponding to list of outliers}
  \item{...}{additional arguments to \code{\link{plot}}}
}
\details{
  Outliers must be determined in advance and stored as a list, with one
  element per chromosome.  The X chromosome must be a list of two
  elements, "M" and "F".  Each element should contain a vector of
  ids corresponding to the row names of \code{mean.intensities}.

  If \code{sep=TRUE}, \code{labels} must also be specified.
  \code{labels} should be a list that corresponds exactly to the elements
  of \code{outliers}.
}
\author{Cathy Laurie}
\seealso{\code{\link{meanIntensityByScanChrom}}}
\examples{
# calculate mean intensity
library(GWASdata)
file <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
nc <- NcdfIntensityReader(file)
data(affy_scan_annot)
scanAnnot <- ScanAnnotationDataFrame(affy_scan_annot)
intenData <- IntensityData(nc, scanAnnot=scanAnnot)
meanInten <- meanIntensityByScanChrom(intenData)
intenMatrix <- meanInten$mean.intensity

# find outliers
outliers <- list()
sex <- scanAnnot$sex
id <- scanAnnot$scanID
allequal(id, rownames(intenMatrix))
for (i in colnames(intenMatrix)) {
  if (i != "X") {
    imean <- intenMatrix[,i]
    imin <- id[imean == min(imean)]
    imax <- id[imean == max(imean)]
    outliers[[i]] <- c(imin, imax)
  } else {
    idf <- id[sex == "F"]
    fmean <- intenMatrix[sex == "F", i]
    fmin <- idf[fmean == min(fmean)]
    fmax <- idf[fmean == max(fmean)]
    outliers[[i]][["F"]] <- c(fmin, fmax)
    idm <- id[sex == "M"]
    mmean <- intenMatrix[sex == "M", i]
    mmin <- idm[mmean == min(mmean)]
    mmax <- idm[mmean == max(mmean)]
    outliers[[i]][["M"]] <- c(mmin, mmax)
  }
}

par(mfrow=c(2,4))
intensityOutliersPlot(intenMatrix, sex, outliers)
}
\keyword{hplot}
