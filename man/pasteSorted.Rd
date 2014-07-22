\name{pasteSorted}
\alias{pasteSorted}
\title{Paste two vectors sorted pairwise}
\description{Read a configuration file}
\usage{
pasteSorted(a, b, sep="/")
}
\arguments{
  \item{a}{vector 1}
  \item{b}{vector 2}
  \item{sep}{a character string to separate the terms.}
}
\value{
  A character vector of the concatenated values, sorted pairwise.
}
\author{Stephanie Gogarten}
\seealso{\code{\link{paste}}}
\examples{
a <- c("A","C","G","T")
b <- c("C","A","T","G")
pasteSorted(a,b)
}
\keyword{manip}
