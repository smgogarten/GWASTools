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
                   chrom.labels = c(1:22,"X","XY","Y","M"),
                   ylim=c(0,1), ylab = "abs(correlation)", ...)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{

  \item{correlations}{A vector of correlations.}
  \item{chromosome}{A vector containing the integer chromosome ID for each SNP.}
  \item{chrom.labels}{A vector of chromosome names to use in the
    plot.}
  \item{ylim}{The limits of the y axis.}
  \item{ylab}{The label for the y axis.}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
  Plots SNP correlations (from, e.g., PCA), versus chromosome.  SNPs are
  evenly spaced along the X axis, so \code{correlations} should
  be in order of position on the chromosome.
}

\author{Cathy Laurie
}


\seealso{\code{\link{manhattanPlot}}
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
correlations <- sample(0.01*(0:100), 100, replace=TRUE)
chromosome <- c(rep(1,50), rep(2,50))
snpCorrelationPlot(correlations, chromosome, chrom.labels=c(1,2))
}


\keyword{hplot }
