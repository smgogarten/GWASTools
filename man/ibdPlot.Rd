\name{ibdPlot}
\alias{ibdPlot}
\alias{ibdAreasDraw}
\alias{ibdAssignRelatedness}
\alias{ibdAssignRelatednessKing}

\title{
Plot theoretical and observed identity by descent values and assign relationships
}

\description{\code{ibdPlot} produces an IBD plot showing observed
  identity by descent values color coded by expected relationship.
  Theoretical boundaries for full-sibling, second-degree, and
  third-degree relatives are plotted in orange.
  \code{ibdAreasDraw} overlays relationship areas for IBD analysis on
  the plot.
  \code{ibdAssignRelatedness} identifies observed relatives.
  \code{ibdAssignRelatedness} identifies observed relatives using the kinship coefficients and IBS0 estimates from the KING model.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usage{
ibdPlot(k0, k1, alpha=0.05, relation=NULL, color=NULL,
        rel.lwd=2, rel.draw=c("FS", "Deg2", "Deg3"), ...)

ibdAreasDraw(alpha=0.05, m=0.04, po.w=0.1, po.h=0.1,
             dup.w=0.1, dup.h=0.1, un.w=0.25, un.h=0.25, rel.lwd=2,
	     xcol=c("cyan","red","blue","lightgreen","magenta","black"))

ibdAssignRelatedness(k0, k1, alpha=0.05, m=0.04, po.w=0.1, po.h=0.1,
                     dup.w=0.1, dup.h=0.1, un.w=0.25, un.h=0.25)

ibdAssignRelatednessKing(ibs0, kc, cut.kc.dup=1/(2^(3/2)),
                         cut.kc.fs=1/(2^(5/2)), cut.kc.deg2=1/(2^(7/2)),
                         cut.kc.deg3=1/(2^(9/2)), cut.ibs0.err=0.003)
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\arguments{
  
  \item{k0}{A vector of k0 values.}
  
  \item{k1}{A vector of k1 values.}
  
  \item{kc}{A vector of kinship coefficient values (KING model).}

  \item{ibs0}{A vector of IBS0 values (KING model).}
  
  \item{alpha}{significance level - finds 100(1-alpha)\% prediction intervals
    for second and third degree relatives and 100(1-alpha)\% prediction
    ellipse for full siblings. 
  }
  \item{relation}{A vector of relationships.  Recognized values are "PO"=parent/offspring, "FS"=full siblings, "HS"=half siblings, "Av"=avuncular, "GpGc"=grandparent-grandchild, "Deg2"=any
    second-degree, "FC"=first
    cousins, "HAv"=half-avuncular, "Deg3"=any third degree,
    "U"=unrelated, and "Q"=unknown.}

  \item{color}{A vector of colors for (k0,k1) points.}

  \item{rel.lwd}{Line width for theoretical full-sib, Deg2, and Deg3 boundaries.}

  \item{rel.draw}{Which theoretical boundaries to plot: one or more of
  "FS" (full-sib), "Deg2" (second-degree), "Deg3" (third-degree).  If
  \code{NULL}, no boundaries are drawn.}
  
  \item{\dots}{Other graphical parameters to pass to \code{\link{plot}} and
  \code{\link{points}}.}

\item{m}{width of rectangle along diagonal line}
\item{po.w}{width of parent-offspring rectangle}
\item{po.h}{height of parent-offspring rectangle}
\item{dup.w}{width of duplicate rectangle}
\item{dup.h}{height of duplicate rectangle}
\item{un.w}{width of unrelated rectangle}
\item{un.h}{height of unrelated rectangle}
\item{xcol}{colors for parent-offspring, full-sib, Deg2, Deg3, dup & unrelated areas}

\item{cut.kc.dup}{Kinship coefficient threshold for dividing duplicates and first degree relatives.}
\item{cut.kc.fs}{Kinship coefficient threshold for dividing full siblings and second degree relatives.}
\item{cut.kc.deg2}{Kinship coefficient threshold for dividing second and third degree relatives.}
\item{cut.kc.deg3}{Kinship coefficient threshold for dividing third degree relatives and unrelated.}
\item{cut.ibs0.err}{IBS0 threshold for dividing parent-offsprings pairs from other relatives. Should be 0, but is usually slightly higher due to genotyping errors.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\details{

\code{ibdPlot} produces an IBD plot showing observed identity by descent
values color coded by expected relationship, typically as deduced from
pedigree data.  Points are plotted according to their corresponding
value in the \code{color} vector, and the \code{relation} vector is used
to make the plot legend.  In addition to the relationships listed above,
any relationships output from \code{\link{pedigreePairwiseRelatedness}}
will be recognized.

Theoretical boundary for full-sibs is indicated by ellipse and
 boundaries for second and third degree intervals are indicated in orange. For full-sibs, 100(1-alpha)\% prediction 
 ellipse is based on assuming bivariate normal distribution with known
 mean and covariance matrix. For second degree (half siblings, avuncular,
 grandparent-grandchild) and third degree (first cousins), 100(1-alpha)\%
 prediction intervals for k1 are based on assuming normal distribution with known mean and
variance, computed as in Hill and Weir (2011).

\code{ibdAreasDraw} overlays relationship areas on the plot to help
with analyzing observed relationships.

\code{ibdAssignRelatedness} identifies relatives based on their (k0, k1) coordinates.

\code{ibdAssignRelatednessKing} identifies relatives based on their (ibs0, kc) coordinates (KING model).

}

\value{
  \code{ibdAssignRelatedness} and \code{ibdAssignRelatednessKing} return a vector of relationships with
  values "Dup"=duplicate, "PO"=parent-offspring, "FS"=full sibling,
  "Deg2"=second degree, "Deg3"=third degree, "U"=unrelated, and
  "Q"=unknown.
  
}

\references{
  Hill, W.G. and B.S. Weir, Variation in actual relationship as a consequence of mendelian sampling and linkage, Genet. Res. Camb. (2011), 93, 47-64.
  
  Manichaikul, A., Mychaleckyj J.C., Rich S.S., Daly K., Sale M., and Chen W.M., Robust relationship inference in genome-wide association studies, Bioinformatics (2010), 26(22), 2867-2873.
}

\author{Cathy Laurie, Cecelia Laurie, and Adrienne Stilp}

\seealso{
  \code{\link{relationsMeanVar}},
  \code{\link{pedigreePairwiseRelatedness}}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
k0 <- c(0, 0, 0.25, 0.5, 0.75, 1)
k1 <- c(0, 1, 0.5, 0.5, 0.25, 0)
exp.rel <- c("Dup", "PO", "FS", "HS", "FC", "U")
ibdPlot(k0, k1, relation=exp.rel)
ibdAreasDraw()
obs.rel <- ibdAssignRelatedness(k0, k1)

kc <- c(0.5, 0.25, 0.25, 0.125, 0.063, 0)
ibs0 <- c(0, 0, 0.25, 0.5, 0.75, 1)
obs.rel.king <- ibdAssignRelatednessKing(ibs0, kc)
}


\keyword{hplot}
\keyword{manip}
