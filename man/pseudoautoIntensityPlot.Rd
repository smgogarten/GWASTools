\name{pseudoautoIntensityPlot}
\alias{pseudoautoIntensityPlot}
\title{ Plot B Allele Frequency and Log R Ratio for the X and Y chromosomes,
  overlaying XY SNPs}



\description{
This function plots X, Y and pseudoautosomal SNPs on BAF/LRR plots.
}



\usage{
pseudoautoIntensityPlot(intenData, scan.ids, code=NULL,
  plotY=FALSE, hg.build=c("hg18", "hg19"), ...)
}



\arguments{
  \item{scan.ids}{ A vector containing the sample indices of the plots.}
  \item{intenData}{\code{\link{IntensityData}} object, must contain
  'BAlleleFreq' and 'LogRRatio'}
  \item{code}{ A character vector containing the titles to be used for
    each plot. If \code{NULL} then the title will be the sample number and
    the chromosome.}
  \item{plotY}{If \code{plotY} is \code{TRUE}, the Y chromosome will be
    plotted in addition to X.}
  \item{hg.build}{Human genome bulid number}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}
}



\details{
  The pseudoautosomal regions are highlighted on the plots
  (PAR1 and PAR2 in gray, XTR in yellow),
  and the X, Y, and XY SNPs are plotted in different colors.  The base
  positions for these regions depend on genome build (\code{hg.build}).
  Currently hg18 and hg19 are supported.
  
  By default the output is a 2-panel plot with LRR and BAF for the X
  chromosome.  if \code{plotY} is \code{TRUE}, the output is a 4-panel
  plot with the Y chromosome plotted as well.
}

\references{
  Ross, Mark. T. et al. (2005), The DNA sequence of the human X chromosome.
  Nature, 434: 325-337.  doi:10.1038/nature03440
}

\author{ Caitlin McHugh  }

\seealso{\code{\link{pseudoautosomal}},
  \code{\link{IntensityData}}, \code{\link{GenotypeData}},
  \code{\link{BAFfromGenotypes}} }

\examples{
library(GWASdata)
data(illumina_scan_annot)
scanAnnot <- ScanAnnotationDataFrame(illumina_scan_annot)
blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
blnc <- NcdfIntensityReader(blfile)
intenData <-  IntensityData(blnc, scanAnnot=scanAnnot)

scanID <- getScanID(scanAnnot, index=1)
pseudoautoIntensityPlot(intenData=intenData, scan.ids=scanID)
close(intenData)
}

\keyword{ hplot }
