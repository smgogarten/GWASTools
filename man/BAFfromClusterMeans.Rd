\name{BAFfromClusterMeans}
\alias{BAFfromClusterMeans}
\title{B Allele Frequency & Log R Ratio Calculation}
\description{
  This function calculates the B allele frequency and the log R ratio
  values from the mean R and theta values for each cluster.  The values are written to a netCDF file which is assumed to exist with proper variables and size.
}


\usage{
BAFfromClusterMeans(intenData, bl.ncdf.filename,
                    clusterMeanVars = c("tAA","tAB","tBB","rAA","rAB","rBB"),
                    verbose = TRUE)
}
\arguments{
  \item{intenData}{\code{\link{IntensityData}} object holding the X and
    Y intensity data from which the B allele frequency and log R ratio
    are calculated. }
  \item{bl.ncdf.filename}{ The filepath for a previously created netCDF file to hold the B allele frequency and log R ratio values. }
  \item{clusterMeanVars}{Character vector indicating the names of the
    cluster mean columns in the SNP annotation of intenData.  Must be in order
    (tAA,tAB,tBB,rAA,rAB,rBB).}
  \item{verbose}{Logical value specifying whether to show progress information.}
}
\details{
  Because this function can take a considerable amount of time and space, sufficient attention should be given to the value used for \code{block.size}.  The file specified by \code{bl.ncdf.filename} is assumed to have variables 'BAlleleFreq' and 'LogRRatio' to which the proper values are written.
}
\value{
  The netCDF file stored in the \code{bl.ncdf.filename} path is populated with values of B allele frequency and the log R ratio at the completion of this function.  
}
\author{ Stephanie Gogarten, Caitlin McHugh }
\references{Peiffer D.A., Le J.M., Steemers F.J., Chang W., Jenniges T., and et al. High-resolution genomic profiling of chromosomal aberrations using infinium whole-genome genotyping. Genome Research, 16:1136-1148, 2006.}
\seealso{\code{\link{IntensityData}}, \code{\link{BAFfromClusterMeans}}}
\examples{
# create IntensityData object from netCDF
library(GWASdata)
xyfile <- system.file("extdata", "illumina_qxy.nc", package="GWASdata")
xyNC <- NcdfIntensityReader(xyfile)
data(illumina_snp_annot)
snpAnnot <- SnpAnnotationDataFrame(illumina_snp_annot)
xyData <- IntensityData(xyNC, snpAnnot=snpAnnot)
nsamp <- nscan(xyData)

# create netCDF file to hold BAF/LRR data
blfile <- tempfile()
ncdfCreate(illumina_snp_annot, blfile, variables=c("BAlleleFreq","LogRRatio"), n.samples=nsamp)

# calculate BAF and LRR
BAFfromClusterMeans(xyData, blfile, verbose=FALSE)

# read output
blNC <- NcdfIntensityReader(blfile)
blData <- IntensityData(blNC)
baf <- getBAlleleFreq(blData)
lrr <- getLogRRatio(blData)

close(xyNC)
close(blNC)
file.remove(blfile)
}
\keyword{ datagen }
\keyword{ manip }
