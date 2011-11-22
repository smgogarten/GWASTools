\name{ibdPlot}
\alias{ibdPlot}
\alias{ibdAreasDraw}
\alias{ibdAssignRelatedness}

\title{
Plot theoretical and observed identity by descent values and assign relationships
}

\description{\code{ibdPlot} produces an IBD plot showing observed
  identity by descent values color coded by expected relationship.
  Theoretical boundaries for full-sib, half-sib, and first-cousins are
 plotted in orange.
  \code{ibdAreasDraw} overlays relationship areas for IBD analysis on
  the plot.
  \code{ibdAssignRelatedness} identifies observed relatives.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usage{
ibdPlot(k0, k1, alpha=0.05, relation=NULL, color=NULL, rel.lwd=2, ...)

ibdAreasDraw(alpha=.05, m=0.04, po.w=0.1, po.h=0.1,
             dup.w=0.1, dup.h=0.1, un.w=0.25, un.h=0.25, rel.lwd=2,
	     xcol=c("cyan","red","blue","lightgreen","magenta","black"))

ibdAssignRelatedness(k0, k1, alpha=0.05, m=0.04, po.w=0.1, po.h=0.1,
                     dup.w=0.1, dup.h=0.1, un.w=0.25, un.h=0.25)
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\arguments{
  
  \item{k0}{A vector of k0 values.}
  
  \item{k1}{A vector of k1 values.}

  \item{alpha}{significance level - finds 100(1-alpha)\% prediction intervals
    for half-sibs and first cousins and 100(1-alpha)\% prediction ellipse for full-sibs 
  }
  \item{relation}{A vector of relationships.}

  \item{color}{A vector of colors for (k0,k1) points.}

  \item{rel.lwd}{Line width for theoretical full-sib, half-sib, and
  first-cousin boundaries.}
  
  \item{\dots}{Other graphical parameters to pass to \code{\link{plot}} and
  \code{\link{points}}.}

\item{m}{width of rectangle along diagonal line}
\item{po.w}{width of parent-offspring rectangle}
\item{po.h}{height of parent-offspring rectangle}
\item{dup.w}{width of duplicate rectangle}
\item{dup.h}{height of duplicate rectangle}
\item{un.w}{width of unrelated rectangle}
\item{un.h}{height of unrelated rectangle}
\item{xcol}{colors for parent-offspring, full-sib, half-sib, first cousin, dup & unrelated areas}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\details{

\code{ibdPlot} produces an IBD plot showing observed identity by descent
values color coded by expected relationship, typically as deduced from
pedigree data.  Points are plotted according to their corresponding
value in the \code{color} vector, and the \code{relation} vector is used
to make the plot legend.

Theoretical boundary for full-sibs is indicated by ellipse and boundaries for half-sib and
first cousin intervals are indicated in orange. For full-sibs, 100(1-alpha)\% prediction 
 ellipse is based on assuming bivariate normal distribution with known
 mean and covariance matrix. For half-sibs and first-cousins, 100(1-alpha)\%
 prediction intervals for k1 are based on assuming normal distribution with known mean and
variance, computed as in Hill and Weir (2011).

\code{ibdAreasDraw} overlays relationship areas on the plot to help
with analyzing observed relationships. For full-sibs, 100(1-alpha)\% prediction 
 ellipse is based on assuming bivariate normal distribution with known
 mean and covariance matrix. For half-sibs and first-cousins, 100(1-alpha)\%
 prediction intervals for k1 are based on assuming normal distribution with known mean and
variance.

\code{ibdAssignRelatedness} identifies relatives based on their (k0, k1) coordinates.

}

\value{
  \code{ibdAssignRelatedness} returns a vector of relationships with
  values "Dup"=duplicate, "PO"=parent-offspring, "FS"=full sibling,
  "HS"=half-sibling-like, "FC"=first cousin, "U"=unrelated, and
  "Q"=unknown.
}

\references{
  Hill, W.G. and B.S. Weir, Variation in actual relationship as a consequence of mendelian sampling and linkage, Genet. Res. Camb. (2011), 93, 47-64.
}

\author{Cathy Laurie and Cecelia Laurie}
\seealso{ \code{\link{relationsMeanVar}} }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
k0 <- c(0, 0, 0.25, 0.5, 0.75, 1)
k1 <- c(0, 1, 0.5, 0.5, 0.25, 0)
exp.rel <- c("Dup", "PO", "FS", "HS", "FC", "U")
ibdPlot(k0, k1, relation=exp.rel)
ibdAreasDraw()
obs.rel <- ibdAssignRelatedness(k0, k1)
}


\keyword{hplot}
\keyword{manip}
