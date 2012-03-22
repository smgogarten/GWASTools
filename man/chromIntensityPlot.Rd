\name{chromIntensityPlot}
\alias{chromIntensityPlot}
\title{ Plot B Allele Frequency and/or Log R Ratio, R or Theta values for samples by probe position on a chromosome}



\description{
This function creates plots for one or more of the 'B AlleleFreq', 'Log R Ratio', 'R' or 'Theta' values for given sample by chromosome combinations.  
}



\usage{
chromIntensityPlot(intenData, scan.ids, chrom.ids, 
                   type = c("BAF/LRR", "BAF", "LRR", "R", "Theta", "R/Theta"), 
                   code = NULL, main.txt = NULL, 
                   abln = NULL,  horizln = c(1/2, 1/3, 2/3), 
                   colorGenotypes = FALSE, genoData = NULL,
                   colorBatch = FALSE, batch.column = NULL,
                   snp.exclude = NULL, cex=0.5, ...)
}



\arguments{
  \item{intenData}{\code{\link{IntensityData}} object, must contain at least one of 'BAlleleFreq', 'LogRRatio', 'X', 'Y'. }
  \item{scan.ids}{ A vector containing the sample indices of the plots.}
  \item{chrom.ids}{ A vector containing the chromosome indices of the plots.}
  \item{type}{ The type of plot to be created.  'BAF/LRR' creates both 'B Allele Freq' and 'Log R Ratio' plots. 'R/Theta' creates both 'R' and 'Theta' plots. }
  \item{code}{ A character vector containing extra information to include in the main title. }
  \item{main.txt}{ Text that will be written in the title on all plots created. If \code{NULL} then the title will be the sample number and the chromosome.}
  \item{abln}{A vector of values that is of length \code{2*length(scan.ids)}.  Each pair of entries specifies where vertical lines will be drawn in each plot.
This is especially useful when drawing the start \& end breakpoints for anomalies, for example. 
Use -1 as one pair value for plots that warrant only one line. By default, no lines will be drawn. }
  \item{horizln}{ A vector containing the y-axis values at which a horizontal line will be drawn in B Allele Frequency plots.}
  \item{colorGenotypes}{A logical value specifying whether to color-code
    the points by called genotype.  if \code{TRUE}, genoData must be given also.}
  \item{genoData}{\code{\link{GenotypeData}} object, required if colorGenotypes=\code{TRUE}.}
  \item{colorBatch}{A logical value specifying whether to color-code the
  points by sample batch (e.g, plate).  If \code{TRUE}, batch.column must also
  be specified.}
  \item{batch.column}{A character string indicating which annotation
	variable in intenData should be used as the batch.}
  \item{snp.exclude}{An integer vector giving the IDs of SNPs to exclude
  from the plot.}
  \item{cex}{cex value for points on the plots}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}
}



\details{
For all plots, a vertical line is drawn every one eigth of the
chromosome. For the Log R Ratio plot, the y-axis has been given the range of (-2,2).
}



\author{ Caitlin McHugh, Cathy Laurie }

\seealso{\code{\link{IntensityData}}, \code{\link{GenotypeData}},
  \code{\link{BAFfromGenotypes}} }

\examples{
library(GWASdata)
data(illuminaScanADF)
blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
blnc <- NcdfIntensityReader(blfile)
intenData <-  IntensityData(blnc, scanAnnot=illuminaScanADF)

genofile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
genonc <- NcdfGenotypeReader(genofile)
genoData <-  GenotypeData(genonc, scanAnnot=illuminaScanADF)

scanID <- getScanID(illuminaScanADF, index=1)
chromIntensityPlot(intenData=intenData, scan.ids=scanID,
                   chrom.ids=22, type="BAF/LRR", code="interesting sample",
                   colorGenotypes=TRUE, genoData=genoData)
close(genoData)
close(intenData)
}

\keyword{ hplot }
