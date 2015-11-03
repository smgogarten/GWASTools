\name{BAFfromGenotypes}
\alias{BAFfromGenotypes}
\title{B Allele Frequency & Log R Ratio Calculation}
\description{
  This function calculates the B allele frequency and the log R ratio values for samples by either plate or by study.
}


\usage{
BAFfromGenotypes(intenData, genoData, 
                 filename, file.type = c("gds", "ncdf"),
                 min.n.genotypes = 2, 
                 call.method = c("by.plate", "by.study"), 
                 plate.name = "plate", 
                 block.size = 5000, 
	         precision="single", compress="ZIP_RA",
                 verbose = TRUE)
}
\arguments{
  \item{intenData}{\code{\link{IntensityData}} object holding the X and Y intensity data from which the B allele frequency and log R ratio are calculated. }
  \item{genoData}{\code{\link{GenotypeData}} object. }
\item{filename}{The name of the genotype GDS or netCDF file to create}

\item{file.type}{The type of file to create ("gds" or "ncdf")}

  \item{min.n.genotypes}{ The minimum number of samples for each genotype at any SNP in order to have non-missing B allele freqency and log R ratio.  Setting this parameter to 2 or a similar value is recommended.}
  \item{call.method}{ If call.method is 'by.plate', the B allele frequency and log R ratio are calculated for samples delineated by plates.  This is the default method.  If call.method is 'by.study', the calculation uses all samples at once.  If a study does not have plate specifications, 'by.study' is the call.method that must be used. }
  \item{plate.name}{ Character string specifying the name of the plate
    variable in intenData or genoData.  By default, the plate.name is simply 'plate' but oftentimes there are variations, such as 'plateID' or 'plate.num'. }
  \item{block.size}{ An integer specifying the number of SNPs to be loaded at one time.  The recommended value is around 1000, but should vary depending on computing power. }
\item{precision}{A character value indicating whether floating point
  numbers should be stored as "double" or "single" precision.}

\item{compress}{The compression level for variables in a GDS file (see
\code{\link{add.gdsn}} for options.}

  \item{verbose}{Logical value specifying whether to show progress information.}
}
\details{
  Because this function can take a considerable amount of time and space, sufficient attention should be given to the value used for \code{block.size}.
}
\author{ Caitlin McHugh }
\references{Peiffer D.A., Le J.M., Steemers F.J., Chang W., Jenniges T., and et al. High-resolution genomic profiling of chromosomal aberrations using infinium whole-genome genotyping. Genome Research, 16:1136-1148, 2006.}
\seealso{\code{\link{IntensityData}}, \code{\link{GenotypeData}},
  \code{\link{chromIntensityPlot}}, \code{\link{BAFfromClusterMeans}}}
\examples{
\dontrun{
# create IntensityData and GenotypeData objects from netCDF
library(GWASdata)
data(affySnpADF)
data(affyScanADF)
nsamp <- nrow(affyScanADF)

xyfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
xyNC <- NcdfIntensityReader(xyfile)
xyData <- IntensityData(xyNC, snpAnnot=affySnpADF, scanAnnot=affyScanADF)

genofile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
genoNC <- NcdfGenotypeReader(genofile)
genoData <- GenotypeData(genoNC, snpAnnot=affySnpADF, scanAnnot=affyScanADF)

# calculate BAF and LRR
blfile <- tempfile()
BAFfromGenotypes(xyData, genoData, blfile, file.type="ncdf", min.n.genotypes=2,
                 call.method="by.plate", plate.name="plate")

blNC <- NcdfIntensityReader(blfile)
baf <- getBAlleleFreq(blNC)
lrr <- getLogRRatio(blNC)

close(xyData)
close(genoData)
close(blNC)
file.remove(blfile)
}
}
\keyword{ datagen }
\keyword{ manip }
