\name{manhattanPlot}
\alias{manhattanPlot}

\title{
Manhattan plot for genome wide association tests
}
\description{
Generates a manhattan plot of the results of a genome wide association test.
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
manhattanPlot(p, chromosome, 
              ylim = NULL, trunc.lines = TRUE,
              signif = 5e-8, thinThreshold=NULL, pointsPerBin=10000, ...)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{

  \item{p}{A vector of p-values.}
  \item{chromosome}{A vector containing the chromosome for each SNP.}
  \item{ylim}{The limits of the y axis.  If NULL, the y axis is \code{(0,
    log10(length(p)) + 4)}.}
  \item{trunc.lines}{Logical value indicating whether to show truncation
    lines.}
  \item{signif}{Genome-wide significance level for plotting horizontal
  line.  If \code{signif=NULL}, no line will be drawn.}
  \item{thinThreshold}{if not \code{NULL},  \code{-log10(pval)} threshold for thinning points.}
  \item{pointsPerBin}{number of points to plot in each bin if \code{thinThreshold} is given. Ignored otherwise.}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
  Plots -log10(p) versus chromosome.  Point size is scaled so that
  smaller p values have larger points.

  \code{p} must have the same length as \code{chromosome} and is
  assumed to be in order of position on each chromosome.  Values within
  each chromosome are evenly spaced along the X axis.

  Plot limits are determined as follows: if \code{ylim} is provided, any
  points with \code{-log10(p) > ylim[2]} are plotted as triangles at the
  maximum y value of the plot.  A line will be drawn to indicate
  trunctation (if \code{trunc.lines == TRUE}, the default).  If
  \code{ylim == NULL}, the maximum y value is defined as
  \code{log10(length(p)) + 4)}.
  
If requested with \code{thinThreshold}, points with \code{-log10(pval) < thinThreshold} are thinned before plotting. All points with \code{-log10(pval) >= thinThreshold} are displayed. P-values with \code{-log10(pval) < thinThreshold} are sampled such that \code{pointsPerBin} points are randomly selected from 10 bins with uniform spacing in \code{-log10(pval)} space.


}

\author{Cathy Laurie, Adrienne Stilp}


\seealso{\code{\link{snpCorrelationPlot}}
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
n <- 1000
pvals <- sample(-log10((1:n)/n), n, replace=TRUE)
chromosome <- c(rep(1,400), rep(2,350), rep("X",200), rep("Y",50))
manhattanPlot(pvals, chromosome, signif=1e-7)
manhattanPlot(pvals, chromosome, thinThreshold=2)
}


\keyword{hplot}
