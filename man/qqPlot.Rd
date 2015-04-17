\name{qqPlot}
\alias{qqPlot}

\title{QQ plot for genome wide assocation studies
}
\description{
Generates a Quantile-Quantile plot for -log10 p-values from genome wide association tests.
}
\usage{
qqPlot(pval, truncate = FALSE, ylim = NULL, thinThreshold = NULL, ci=TRUE, ...)
}

\arguments{
  \item{pval}{Vector of p-values
}
  \item{truncate}{Either a logical value indicating whether the y-axis should be truncted to the same range as the x-axis, or a numeric value indicating where to truncate the y-axis. See details.
  }
  \item{ylim}{Limits for the y axis.  Ignored if \code{truncate=TRUE} or \code{truncate} is numeric.}
  \item{thinThreshold}{if not \code{NULL},  \code{-log10(pval)} threshold for thinning points.}
  \item{ci}{logical indicator for whether to add confidence intervals to plots?}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}

}


\details{
The function generates a Quantile-Quantile plot of p-values on a -log10
scale, with the option of truncating the y-axis to the range of the
x-axis \code{(0, -log10(1/length(pval))}.  If the y-axis is truncated,
then points off the top of the plot are denoted by triangles at the upper edge.
The 95\% confidence interval is shaded in gray.

If truncate is set to a numeric value, then \code{ylim} is set to \code{c(0, truncate)} only if the value of truncate is bigger than the maximum \code{-log10(pval)}. (Use the \code{ylim} argument if alternatve behavior is desired.)

If requested with \code{thinThreshold}, points with p-values \code{< -log10(thinThreshold)} are thinned before plotting. All points with \code{-log10(pval) >= thinThreshold} plus 10,000 points with \code{-log10(pval) < thinThreshold} (randomly selected in uniformly-spaced bins of -log10(pval)) are displayed.
}
\author{Cathy Laurie, Matthew P. Conomos, Adrienne Stilp
}

\examples{
pvals <- seq(0, 1, 0.001)
qqPlot(pvals)
qqPlot(pvals, thinThreshold=2)
qqPlot(pvals, truncate=TRUE)
qqPlot(pvals, truncate=10)
}
\keyword{ hplot }
