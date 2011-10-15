\name{qqPlot}
\alias{qqPlot}

\title{QQ plot for genome wide assocation studies
}
\description{
Generates a Quantile-Quantile plot for -log10 p-values from genome wide association tests.
}
\usage{
qqPlot(pval, truncate = FALSE, sub = NULL, ...)
}

\arguments{
  \item{pval}{Vector of p-values
}
  \item{truncate}{A logical value indicating whether the y-axis should be truncted to the same range as the x-axis.
  }
  \item{sub}{A character string to print under the x-axis.}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}
  
}


\details{
The function generates a Quantile-Quantile plot of p-values on a -log10
scale, with the option of truncating the y-axis to the range of the
x-axis \code{(0, -log10(1/length(pval))}.  If the y-axis is truncated, 
then points off the top of the plot are denoted by triangles at the upper edge.  
The 95\% confidence interval is shaded in gray.

If \code{sub=NULL} (the default), the genomic inflation factor lambda is
  calculated from the p-values and printed.
}
\author{Cathy Laurie, Matt Conomos
}

\examples{
pvals <- seq(0, 1, 0.001)
qqPlot(pvals)
}
\keyword{ hplot }
