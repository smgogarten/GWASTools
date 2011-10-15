\name{centromeres}
\alias{centromeres}
\alias{centromeres.hg18}
\alias{centromeres.hg19}
\docType{data}
\title{
Centromere base positions
}
\description{
Centromere base positions from the GRCh36/hg18 and GRCh37/hg19 genome builds.
}
\usage{
data(centromeres.hg18)
data(centromeres.hg19)
}
\format{
  A data frame with the following columns.
  \describe{
    \item{\code{chrom}}{chromosome (1-22, X, Y)}
    \item{\code{left.base}}{starting base position of centromere}
    \item{\code{right.base}}{ending base position of centromere}
  }
}
\source{
UCSC genome browser (\url{http://genome.ucsc.edu}).
}
\note{
  The UCSC genome browser lists two regions for the Y chromosome
  centromere in build hg18.  We removed the positions (12208578,
  12308578) from the centromere table to avoid problems with duplicate
  entries in the code.
}
\examples{
data(centromeres.hg18)
data(centromeres.hg19)
}
\keyword{datasets}
