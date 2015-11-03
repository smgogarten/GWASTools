\name{BAFfromClusterMeans}
\alias{BAFfromClusterMeans}
\title{B Allele Frequency & Log R Ratio Calculation}
\description{
  This function calculates the B allele frequency and the log R ratio
  values from the mean R and theta values for each cluster.
}


\usage{
BAFfromClusterMeans(intenData, filename, file.type = c("gds", "ncdf"),
                    clusterMeanVars = c("tAA","tAB","tBB","rAA","rAB","rBB"),
	            precision="single", compress="ZIP_RA",
                    verbose = TRUE)
}
\arguments{
  \item{intenData}{\code{\link{IntensityData}} object holding the X and
    Y intensity data from which the B allele frequency and log R ratio
    are calculated. }
\item{filename}{The name of the genotype GDS or netCDF file to create}

\item{file.type}{The type of file to create ("gds" or "ncdf")}

  \item{clusterMeanVars}{Character vector indicating the names of the
    cluster mean columns in the SNP annotation of intenData.  Must be in order
    (tAA,tAB,tBB,rAA,rAB,rBB).}

\item{precision}{A character value indicating whether floating point
  numbers should be stored as "double" or "single" precision.}

\item{compress}{The compression level for variables in a GDS file (see
\code{\link{add.gdsn}} for options.}

  \item{verbose}{Logical value specifying whether to show progress information.}
}
\details{
  This function calculates the B allele frequency and the log R ratio
  values from the mean R and theta values for each cluster and writes them to a GDS or NetCDF file.
}
\author{ Stephanie Gogarten, Caitlin McHugh }
\references{Peiffer D.A., Le J.M., Steemers F.J., Chang W., Jenniges T., and et al. High-resolution genomic profiling of chromosomal aberrations using infinium whole-genome genotyping. Genome Research, 16:1136-1148, 2006.}
\seealso{\code{\link{IntensityData}}, \code{\link{BAFfromClusterMeans}}}
\examples{
# create IntensityData object from GDS
library(GWASdata)
xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
xy <- GdsIntensityReader(xyfile)
data(illuminaSnpADF)
xyData <- IntensityData(xy, snpAnnot=illuminaSnpADF)

# calculate BAF and LRR and store in GDS file
blfile <- tempfile()
BAFfromClusterMeans(xyData, blfile, file.type="gds", verbose=FALSE)

# read output
bl <- GdsIntensityReader(blfile)
baf <- getBAlleleFreq(bl)
lrr <- getLogRRatio(bl)

close(xy)
close(bl)
file.remove(blfile)
}
\keyword{ datagen }
\keyword{ manip }
