\name{allequal}
\alias{allequal}
\title{Test if two objects have the same elements}
\description{\code{allequal} tests if two objects have all the same
  elements, including whether they have \code{NA}s in the same place.
}
\usage{allequal(x, y)}
\arguments{
  \item{x}{first object to compare}
  \item{y}{second object to compare}
}
\details{  Unlike
  \code{all(x == y)}, \code{allequal} will return \code{FALSE} if either object is \code{NULL}.  Does not
  check class types, so \code{allequal} will return \code{TRUE} in some cases
  where \code{identical} will return \code{FALSE} (e.g. if two objects are
  identical when coerced to the same class).  \code{allequal} always
  retuns a logical value, so it can be used safely in \code{if} expressions.
}
\value{
  Returns \code{TRUE} if x and y exist and all elements are equal, \code{FALSE} if
  some elements are unequal.
  If there are \code{NA} values, returns \code{TRUE} if is.na(x) == is.na(y) and
  all other elements are equal.
  Returns \code{FALSE} if is.na(x) != is.na(y).
  Retuns \code{FALSE} if x or y (but not both) is \code{NULL}.
}
\author{Stephanie Gogarten}
\seealso{\code{\link{identical}}, \code{\link{all}}, \code{\link{all.equal}}
}
\examples{
x <- c(1,2,NA,4); y <- c(1,2,NA,4);
allequal(x, y) ## TRUE
allequal(1, as.integer(1)) ## TRUE
allequal(1, "1") ## TRUE
}
\keyword{logic}
