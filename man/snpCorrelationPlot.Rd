\name{snpCorrelationPlot}
\alias{snpCorrelationPlot}

\title{
SNP correlation plot
}
\description{
Plots SNP correlation versus chromosome.
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
snpCorrelationPlot(correlations, chromosome, 
                   ylim=c(0,1), ylab = "abs(correlation)", ...)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{

  \item{correlations}{A vector of correlations.}
  \item{chromosome}{A vector containing the chromosome for each SNP.}
  \item{ylim}{The limits of the y axis.}
  \item{ylab}{The label for the y axis.}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
  Plots SNP correlations (from, e.g., PCA), versus chromosome.  
  
  \code{correlations} must have the same length as \code{chromosome} and is
  assumed to be in order of position on each chromosome.  Values within
  each chromosome are evenly spaced along the X axis.

}

\author{Cathy Laurie
}


\seealso{\code{\link{manhattanPlot}}
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
correlations <- sample(0.01*(0:100), 100, replace=TRUE)
chromosome <- c(rep(1,25), rep(2,25), rep("X",25), rep("Y",25))
snpCorrelationPlot(correlations, chromosome)
}


\keyword{hplot }
