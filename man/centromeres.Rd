\name{centromeres}
\alias{centromeres}
\alias{centromeres.hg18}
\alias{centromeres.hg19}
\alias{centromeres.hg38}
\docType{data}
\title{
Centromere base positions
}
\description{
Centromere base positions from the GRCh36/hg18, GRCh37/hg19 and GRCh38/hg38 genome builds.
}
\usage{
data(centromeres.hg18)
data(centromeres.hg19)
data(centromeres.hg38)
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
hg18 and hg19: UCSC genome browser (\url{http://genome.ucsc.edu})

hg38: Genome Reference Consortium (\url{http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/}).
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
data(centromeres.hg38)
}
\keyword{datasets}
