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
              chrom.labels = c(1:22,"X","XY","Y","M"),
              ylim = NULL, trunc.lines = TRUE,
              signif = 5e-8, ...)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{

  \item{p}{A vector of p-values.}
  \item{chromosome}{A vector containing the integer chromosome ID for each SNP.}
  \item{chrom.labels}{A vector of chromosome names to use in the
    plot.}
  \item{ylim}{The limits of the y axis.  If NULL, the y axis is \code{(0,
    log10(length(p)) + 4)}.}
  \item{trunc.lines}{Logical value indicating whether to show truncation
    lines.}
  \item{signif}{Genome-wide significance level for plotting horizontal
  line.  If \code{signif=NULL}, no line will be drawn.}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
  Plots -log10(p) versus chromosome.  Point size is scaled so that
  smaller p values have larger points.

  Plot limits are determined as follows: if \code{ylim} is provided, any
  points with \code{-log10(p) > ylim[2]} are plotted as triangles at the
  maximum y value of the plot.  A line will be drawn to indicate
  trunctation (if \code{trunc.lines == TRUE}, the default).  If
  \code{ylim == NULL}, the maximum y value is defined as
  \code{log10(length(p)) + 4)}.
}

\author{Cathy Laurie
}


\seealso{\code{\link{snpCorrelationPlot}}
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
n <- 1000
pvals <- sample(-log10((1:n)/n), n, replace=TRUE)
chromosome <- c(rep(1,500), rep(2,500))
manhattanPlot(pvals, chromosome, chrom.labels=c(1,2), signif=1e-7)
}


\keyword{hplot}
