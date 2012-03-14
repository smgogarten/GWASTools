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
    \item{\code{chrom}}{chromosome (X or Y)}2
    \item{\code{region}}{region (XTR, PAR1, or PAR2)}
    \item{\code{start.base}}{starting base position of region}
    \item{\code{end.base}}{ending base position of region}
  }
}
\details{
  The XTR region on X is defined as DXS1217 to DXS3.
  The XTR region on Y is defined as SY20 to DXYS1.
}
\source{
UCSC genome browser (\url{http://genome.ucsc.edu}).
}
\references{
Ross, Mark. T. et al. (2005), The DNA sequence of the human X chromosome.
Nature, 434: 325-337.  doi:10.1038/nature03440

Mumm, S., Molini, B., Terrell, J., Srivastava, A., and Schlessinger, D.
(1997), Evolutionary features of the 4-Mb Xq21.3 XY homology region revealed
by a map at 60-kb resolution. Genome Res. 7: 307-314.
}
\examples{
data(pseudoautosomal.hg18)
data(pseudoautosomal.hg19)
}
\keyword{datasets}
