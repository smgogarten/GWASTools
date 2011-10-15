\name{pseudoautosomal}
\alias{pseudoautosomal}
\alias{pseudoautosomal.hg18}
\alias{pseudoautosomal.hg19}
\docType{data}
\title{Pseudoautosomal region base positions}
\description{
  Pseudoautosomal region (XTR, PAR1, PAR2) base positions for the X and
  Y chromsosomes from the GRCh36/hg18 and GRCh37/hg19 genome builds.
}
\usage{
pseudoautosomal.hg18
pseudoautosomal.hg19
}
\format{
  A data.frame with the following columns.
  \describe{
    \item{\code{chrom}}{chromsome (X or Y)}
    \item{\code{region}}{region (XTR, PAR1, or PAR2)}
    \item{\code{start.base}}{starting base position of region}
    \item{\code{end.base}}{ending base position of region}
  }
}
\source{
  UCSC genome browser (\url{http://genome.ucsc.edu}).
}
\references{
  Ross, Mark. T. et al. (2005), The DNA sequence of the human X chromosome.
  Nature, 434: 325-337.  doi:10.1038/nature03440
}
\examples{
data(pseudoautosomal.hg18)
data(pseudoautosomal.hg19)
}
\keyword{datasets}
