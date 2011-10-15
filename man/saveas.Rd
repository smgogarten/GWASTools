\name{saveas}
\alias{saveas}
\title{Save an R object with a new name}
\description{Saves an R object as name in an Rdata file called path/name.RData.}
\usage{saveas(obj, name, path=".")}
\arguments{
  \item{obj}{R object to save}
  \item{name}{character string with the new name for the R object}
  \item{path}{path for the Rdata file (saved file will be path/name.RData)}
}
\details{The suffix ".RData" will be appended to the new object name to
  create the file name, and the file will be written to the \code{path}
  directory.}
\author{Stephanie Gogarten}
\seealso{\code{\link{getobj}}}
\examples{
x <- 1:10
path <- tempdir()
saveas(x, "myx", path)
newfile <- paste(path, "/myx", ".RData", sep="")
load(newfile) # myx now loaded
unlink(newfile)
}
\keyword{utilities}
