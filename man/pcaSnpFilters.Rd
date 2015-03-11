\name{pcaSnpFilters}
\alias{pcaSnpFilters}
\alias{pcaSnpFilters.hg18}
\alias{pcaSnpFilters.hg19}
\alias{pcaSnpFilters.hg38}
\docType{data}
\title{Regions of SNP-PC correlation to filter for Principal Component Analysis}
\description{
  Base positions for the LCT (2q21), HLA (including MHC), and inversion (8p23, 17q21.31)
  regions from the GRCh36/hg18, GRCh37/hg19 and GRCh38/hg38 genome genome builds.
}
\usage{
pcaSnpFilters.hg18
pcaSnpFilters.hg19
pcaSnpFilters.hg38
}
\format{
  A data.frame with the following columns.
  \describe{
    \item{\code{chrom}}{chromsome}
    \item{\code{start.base}}{starting base position of region}
    \item{\code{end.base}}{ending base position of region}
    \item{\code{comment}}{description of the region}
  }
}
\details{
These regions result in high SNP-PC
  correlation if they are included in Principal Component Analysis
  (PCA).  The \code{pcaSnpFilters} datasets can be used to filter SNPs prior to running PCA
  to avoid correlations.
}
\source{
  UCSC genome browser (\url{http://genome.ucsc.edu}).
}
\references{
  Novembre, John et al. (2008), Genes mirror geography within Europe.
  Nature, 456: 98-101.  doi:10.1038/nature07331
}
\seealso{
  \code{\link{snpCorrelationPlot}}, \code{SNPRelate}
}
\examples{
data(pcaSnpFilters.hg18)
data(pcaSnpFilters.hg19)
data(pcaSnpFilters.hg38)
}
\keyword{datasets}
