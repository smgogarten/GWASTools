\name{readWriteFirst}
\alias{readWriteFirst}
\title{Read and write the first n lines of a file}

\description{
  Read first \code{n} lines of \code{filein} and write them
  to \code{fileout}, where \code{filein} and \code{fileout} are file
  names.
}

\usage{
readWriteFirst(filein, fileout, n)
}

\arguments{
  \item{filein}{input file}
  \item{fileout}{output file}
  \item{n}{number of lines to write}
}

\author{Cathy Laurie}

\examples{
path <- system.file("extdata", "affy_raw_data", package="GWASdata")
file <- paste(path, list.files(path)[1], sep="/")
outf <- tempfile()
readWriteFirst(file, outf, 20)
file.remove(outf)
}

\keyword{IO}
\keyword{file}
