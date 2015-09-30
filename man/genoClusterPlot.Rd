\name{genoClusterPlot}
\alias{genoClusterPlot}
\alias{genoClusterPlotByBatch}
\title{SNP cluster plots}
\description{
  Generates either X,Y or R,Theta cluster plots for specified SNP's.
}


\usage{
genoClusterPlot(intenData, genoData, plot.type = c("RTheta", "XY"), 
                snpID, main.txt = NULL, by.sex = FALSE,
                scan.sel = NULL, scan.hilite = NULL,
                start.axis.at.0 = FALSE,
		colors = c("default", "neon", "primary"),
                verbose = TRUE, ...)

genoClusterPlotByBatch(intenData, genoData, plot.type = c("RTheta", "XY"), 
                       snpID, batchVar, main.txt = NULL, scan.sel = NULL, 
                       colors = c("default", "neon", "primary"),
                       verbose = TRUE, ...)
}
\arguments{
  \item{intenData}{\code{\link{IntensityData}} object containing 'X' and 'Y' values.}
  \item{genoData}{\code{\link{GenotypeData}} object}
  \item{plot.type}{The type of plots to generate.  Possible values are
                       "RTheta" (default) or "XY".}
  \item{snpID}{A numerical vector containing the SNP number for each plot. }
  \item{batchVar}{A character string indicating which annotation
	variable should be used as the batch.}
  \item{main.txt}{A character vector containing the title to give to each plot.}
  \item{by.sex}{Logical value specifying whether to indicate sex on the
	plot.  If \code{TRUE}, sex must be present in intenData or genoData.}
  \item{scan.sel}{integer vector of scans to include in the plot.  If \code{NULL},
	all scans will be included.}
  \item{scan.hilite}{integer vector of scans to highlight in the plot
	with different colors.  If \code{NULL}, all scans will be plotted with
	the same colors.}
  \item{start.axis.at.0}{Logical for whether the min value of each axis
                       should be 0.}
  \item{colors}{Color scheme to use for genotypes. "default" is colorblind safe (colorbrewer Set2), "neon" is bright orange/green/fuschia, and "primary" is red/green/blue.}
  \item{verbose}{Logical value specifying whether to show progress.}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}
}
\details{
  Either 'RTheta' (default) or 'XY' plots can be generated.  R and Theta
  values are computed from X and Y using the formulas \code{r <- x+y} and
  \code{theta <- atan(y/x)*(2/pi)}.
      
  If \code{by.sex==TRUE}, females are indicated with circles and males
  with crosses.
}
\author{ Caitlin McHugh }
\seealso{\code{\link{IntensityData}}, \code{\link{GenotypeData}}}
\examples{
# create data object
library(GWASdata)
data(illuminaScanADF, illuminaSnpADF)

xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
xy <- GdsIntensityReader(xyfile)
xyData <-  IntensityData(xy, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
geno <- GdsGenotypeReader(genofile)
genoData <-  GenotypeData(geno, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

# select first 9 snps
snpID <- illuminaSnpADF$snpID[1:9]
rsID <- illuminaSnpADF$rsID[1:9]

par(mfrow=c(3,3)) # plot 3x3
genoClusterPlot(xyData, genoData, snpID=snpID, main.txt=rsID)

# select samples
scan.sel <- illuminaScanADF$scanID[illuminaScanADF$race == "CEU"]
genoClusterPlot(xyData, genoData, snpID=snpID, main.txt=rsID,
                scan.sel=scan.sel, by.sex=TRUE)

genoClusterPlot(xyData, genoData, snpID=snpID, main.txt=rsID,
                scan.hilite=scan.sel)
close(xyData)
close(genoData)

## affy data - cluster plots by plate
data(affyScanADF, affySnpADF)

xyfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
xy <- NcdfIntensityReader(xyfile)
xyData <-  IntensityData(xy, scanAnnot=affyScanADF, snpAnnot=affySnpADF)

genofile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
geno <- NcdfGenotypeReader(genofile)
genoData <-  GenotypeData(geno, scanAnnot=affyScanADF, snpAnnot=affySnpADF)

# select first 9 snps
snpID <- affySnpADF$snpID[1:9]
rsID <- affySnpADF$rsID[1:9]

genoClusterPlotByBatch(xyData, genoData, snpID=snpID, main.txt=rsID,
                       batchVar="plate")
close(xyData)
close(genoData)

}
\keyword{ hplot }
