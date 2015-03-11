\name{HLA}
\alias{HLA}
\alias{HLA.hg18}
\alias{HLA.hg19}
\alias{HLA.hg38}
\docType{data}
\title{HLA region base positions}
\description{
  HLA region base positions from the GRCh36/hg18, GRCh37/hg19 and GRCh38/hg38 genome builds.
}
\usage{
HLA.hg18
HLA.hg19
HLA.hg38
}
\format{
  A data.frame with the following columns.
  \describe{
    \item{\code{chrom}}{chromsome}
    \item{\code{start.base}}{starting base position of region}
    \item{\code{end.base}}{ending base position of region}
  }
}
\source{
  UCSC genome browser (\url{http://genome.ucsc.edu}).
}
\references{
  Mehra, Narinder K. and Kaur, Gurvinder (2003), MHC-based vaccination approaches: progress and perspectives.
  Expert Reviews in Molecular Medicine, Vol. 5: 24.  doi:10.1017/S1462399403005957
}
\examples{
data(HLA.hg18)
data(HLA.hg19)
data(HLA.hg38)
}
\keyword{datasets}
