\name{anomSegStats}
\alias{anomSegStats}
\alias{anomStatsPlot}
\title{Calculate LRR and BAF statistics for anomalous segments}

\description{
  Calculate LRR and BAF statistics for anomalous segments and plot results
}

\usage{
anomSegStats(intenData, genoData, snp.ids, anom, centromere,
  lrr.cut = -2, verbose = TRUE)

anomStatsPlot(intenData, genoData, anom.stats, snp.ineligible,
  plot.ineligible = FALSE, centromere = NULL,
  brackets = c("none", "bases", "markers"), brkpt.pct = 10,
  whole.chrom = FALSE, win = 5, win.calc = FALSE, win.fixed = 1,
  zoom = c("both", "left", "right"), main = NULL, info = NULL,
  ideogram = TRUE, ideo.zoom = FALSE, ideo.rect = TRUE,
  mult.anom = FALSE, cex = 0.5, cex.leg = 1.5, ...)
}

\arguments{
  \item{intenData}{
  An \code{\link{IntensityData}} object containing BAlleleFreq and LogRRatio.
 The order of the rows of intenData and the snp annotation
  are expected to be by chromosome and then by position within chromosome.
}
  \item{genoData}{
  A \code{\link{GenotypeData}} object.  The order of the rows of intenData and the snp annotation
  are expected to be by chromosome and then by position within chromosome.  
}
  \item{snp.ids}{
  vector of eligible SNP ids.  Usually exclude failed and intensity-only SNPS.
  Also recommended to exclude an HLA region on chromosome 6 and
  XTR region on X chromosome.  See \code{\link{HLA}} and \code{\link{pseudoautosomal}}.
  If there are SNPs annotated in the centromere gap, exclude these as
  well (see \code{\link{centromeres}}).
x}
\item{anom}{data.frame of detected chromosome anomalies.  Names must include "scanID",
  "chromosome", "left.index", "right.index", "sex", "method", "anom.id".
  Valid values for "method" are "BAF" or "LOH" referring to whether the anomaly
was detected by BAF method (\code{\link{anomDetectBAF}}) or by LOH method
 (\code{\link{anomDetectLOH}}).   
  Here "left.index" and "right.index" are row indices of intenData with left.index < right.index.
}
\item{centromere}{
  data.frame with centromere position info.  Names must include
  "chrom", "left.base", "right.base".  Valid values for "chrom" are
  1:22, "X", "Y", "XY".  Here "left.base" and "right.base"
  are start and end base positions of the centromere location,
  respectively.  Centromere data tables are provided in \code{\link{centromeres}}.
}
\item{lrr.cut}{
  count the number of eligible LRR values less than \code{lrr.cut}
}
\item{verbose}{
  whether to print the scan id currently being processed
}
\item{anom.stats}{
  data.frame of chromosome anomalies with statistics, usually the output
  of \code{anomSegStats}.  Names must include "anom.id", "scanID", "chromosome",
  "left.index", "right.index", "method", "nmark.all", "nmark.elig", "left.base", "right.base",
  "nbase", "non.anom.baf.med", "non.anom.lrr.med", "anom.baf.dev.med",
  "anom.baf.dev.5", "anom.lrr.med", "nmark.baf", "nmark.lrr". Left and right refer
 to start and end, respectively, of the anomaly, in position order.
}
  \item{snp.ineligible}{
  vector of ineligible snp ids (e.g., intensity-only, failed snps, XTR and HLA regions).
 See \code{\link{HLA}} and \code{\link{pseudoautosomal}}.
}
\item{plot.ineligible}{
  whether or not to include ineligible points in the plot for LogRRatio
}
\item{brackets}{
  type of brackets to plot around breakpoints - none, use base length, use number of markers (note that using markers give asymmetric brackets);
 could be used, along with \code{brkpt.pct}, to assess general accuracy of end points of the anomaly
}
\item{brkpt.pct}{
  percent of anomaly length in bases (or number of markers) for width of brackets
}
\item{whole.chrom}{
  logical to plot the whole chromosome or not (overrides \code{win} and \code{zoom})
}
\item{win}{
  size of the window (a multiple of anomaly length) surrounding the anomaly to plot
}
\item{win.calc}{
  logical to calculate window size from anomaly length; overrides \code{win} and gives window of fixed length given by \code{win.fixed}
}
\item{win.fixed}{
  number of megabases for window size when \code{win.calc=TRUE}
}
\item{zoom}{
  indicates whether plot includes the whole anomaly ("both") or zooms on just the left or right breakpoint; "both" is default
}
\item{main}{ Vector of titles for upper (LRR) plots. If \code{NULL}, titles will
  include anom.id, scanID, sex, chromosome, and detection method. }
\item{info}{
  character vector of extra information to include in the main title of
  the upper (LRR) plot
}
\item{ideogram}{logical for whether to plot a chromosome ideogram under
  the BAF and LRR plots.}
  \item{ideo.zoom}{logical for whether to zoom in on the ideogram to
    match the range of the BAF/LRR plots}
  \item{ideo.rect}{logical for whether to draw a rectangle on the
    ideogram indicating the range of the BAF/LRR plots}
  \item{mult.anom}{logical for whether to plot multiple anomalies from
  the same scan-chromosome pair on a single plot.  If \code{FALSE}
  (default), each anomaly is shown on a separate plot.}
\item{cex}{
  cex value for points on the plots
}
  \item{cex.leg}{cex value for the ideogram legend}
  \item{\dots}{Other parameters to be passed directly to \code{\link{plot}}.}
}

\details{ \code{anomSegStats} computes various statistics of the input anomalies.
 Some of these are basic statistics for the characteristics of the anomaly and for measuring deviation of LRR or BAF from expected.
 Other statistics are used in downstrean quality control analysis, including detecting
 terminal anomalies and investigating centromere-spanning anomalies. 


\code{anomStatsPlot} produces separate png images of each anomaly in \code{anom.stats}. Each image consists of
an upper plot of LogRRatio values and a lower plot of BAlleleFrequency values for 
a zoomed region around the anomaly or whole chromosome (depending up parameter
choices). Each plot has vertical lines demarcating the anomaly and horizontal lines
displaying certain statistics from \code{anomSegStats}. The upper plot title
includes sample number and chromosome. Further plot annotation describes which
anomaly statistics are represented.

}

\value{
  \code{anomSegStats} produces a data.frame with the variables for \code{anom} plus the following columns: 
Left and right refer to position order with left < right.


\item{nmark.all}{total number of SNP markers on the array from left.index to right.index inclusive}
\item{nmark.elig}{total number of eligible SNP markers on the array from left.index to right.index, inclusive.  
See \code{snp.ids} for definition of eligible SNP markers.}
\item{left.base}{base position corresponding to left.index}
\item{right.base}{base position corresponding to right.index}
\item{nbase}{number of bases from left.index to right.index, inclusive}
\item{non.anom.baf.med}{BAF median of non-anomalous segments on all autosomes for the associated sample,
  using eligible heterozygous or missing SNP markers }
\item{non.anom.lrr.med}{LRR median of non-anomalous segments on all autosomes for the associated sample,
  using eligible SNP markers }
\item{non.anom.lrr.mad}{MAD for LRR of non-anomalous segments on all autosomes for the associated sample, using eligible SNP markers}
\item{anom.baf.dev.med}{BAF median of deviations from \code{non.anom.baf.med} of points used to detect anomaly (eligible and heterozygous or missing) 
 }
\item{anom.baf.dev.5}{median of BAF deviations from 0.5, using eligible heterozygous or missing SNP markers in anomaly }
\item{anom.baf.dev.mean}{mean of BAF deviations from \code{non.anom.baf.med}, using eligible heterozygous or missing SNP markers in anomaly }
\item{anom.baf.sd}{standard deviation of BAF deviations from \code{non.anom.baf.med}, using eligible heterozygous or missing SNP markers in anomaly}
\item{anom.baf.mad}{MAD of BAF deviations from \code{non.anom.baf.med}, using eligible heterozygous or missing SNP markers in anomaly }

\item{anom.lrr.med}{LRR median of eligible SNP markers within the anomaly}
\item{anom.lrr.sd}{standard deviation of LRR for eligible SNP markers within the anomaly}
\item{anom.lrr.mad}{MAD of LRR for eligible SNP markers within the anomaly}
\item{nmark.baf}{number of SNP markers within the anomaly eligible for BAF detection (eligible markers that are heterozygous or missing)}
\item{nmark.lrr}{number of SNP markers within the anomaly eligible for LOH detection (eligible markers)}

\item{cent.rel}{position relative to centromere - left, right, span}

\item{left.most}{T/F for whether the anomaly is the left-most anomaly for this sample-chromosome,
  i.e. no other anomalies with smaller start base position }
\item{right.most}{T/F whether the anomaly is the right-most anomaly for this sample-chromosome,
  i.e. no other anomalies with larger end base position}
\item{left.last.elig}{T/F for whether the anomaly contains the last eligible SNP marker going to the left (decreasing position)
}
\item{right.last.elig}{T/F for whether the anomaly contains the last eligible SNP marker going to the right (increasing position)
}

\item{left.term.lrr.med}{median of LRR for all eligible SNP markers from left-most eligible marker to the left telomere
 (only calculated for the most distal anom)}
\item{right.term.lrr.med}{median of LRR for all eligible markers from right-most eligible marker to the right telomere
 (only calculated for the most distal anom)}
\item{left.term.lrr.n}{sample size for calculating \code{left.term.lrr.med}}
\item{right.term.lrr.n}{sample size for calculating \code{right.term.lrr.med}}

\item{cent.span.left.elig.n}{number of eligible markers on the left side of centromere-spanning anomalies
}
\item{cent.span.right.elig.n}{number of eligible markers on the right side of centromere-spanning anomalies
}
\item{cent.span.left.bases}{length of anomaly (in bases) covered by eligible markers on the left side of the centromere}
\item{cent.span.right.bases}{length of anomaly (in bases) covered by eligible markers on the right side of the centromere}
\item{cent.span.left.index}{index of eligible marker left-adjacent to centromere;
 recall that index refers to row indices of \code{intenData}}
\item{cent.span.right.index}{index of elig marker right-adjacent to centromere}

\item{bafmetric.anom.mean}{mean of BAF-metric values within anomaly, using eligible heterozygous or missing SNP markers   BAF-metric values were used in the
  detection of anomalies. See \code{\link{anomDetectBAF}} for definition of BAF-metric}

\item{bafmetric.non.anom.mean}{mean of BAF-metric values within non-anomalous segments 
  across all autosomes for the associated sample, using eligible heterozygous or missing SNP markers } 

\item{bafmetric.non.anom.sd}{standard deviation of BAF-metric values within non-anomalous segments 
  across all autosomes for the associated sample, using eligible heterozygous or missing SNP markers }

\item{nmark.lrr.low}{number of eligible markers within anomaly with LRR values less than \code{lrr.cut}}
}

\author{Cathy Laurie}
\note{The non-anomalous statistics are computed over all autosomes for 
 the sample associated with an anomaly. Therefore the accuracy of these statistics
 relies on the input anomaly data.frame including all autosomal anomalies for a given sample.

}

\seealso{
  \code{\link{anomDetectBAF}}, \code{\link{anomDetectLOH}}
}

\examples{
library(GWASdata)
data(illuminaScanADF, illuminaSnpADF)

blfile <- system.file("extdata", "illumina_bl.gds", package="GWASdata")
bl <- GdsIntensityReader(blfile)
blData <-  IntensityData(bl, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
geno <- GdsGenotypeReader(genofile)
genoData <-  GenotypeData(geno, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

scan.ids <- illuminaScanADF$scanID[1:2]
chrom.ids <- unique(illuminaSnpADF$chromosome)
snp.ids <- illuminaSnpADF$snpID[illuminaSnpADF$missing.n1 < 1]
snp.failed <- illuminaSnpADF$snpID[illuminaSnpADF$missing.n1 == 1]

# example results from anomDetectBAF
baf.anoms <- data.frame("scanID"=rep(scan.ids[1],2), "chromosome"=rep(21,2),
  "left.index"=c(100,300), "right.index"=c(200,400), sex=rep("M",2),
  method=rep("BAF",2), anom.id=1:2, stringsAsFactors=FALSE)

# example results from anomDetectLOH
loh.anoms <- data.frame("scanID"=scan.ids[2],"chromosome"=22,
  "left.index"=400,"right.index"=500, sex="F", method="LOH",
  anom.id=3, stringsAsFactors=FALSE)

anoms <- rbind(baf.anoms, loh.anoms)
data(centromeres.hg18)
stats <- anomSegStats(blData, genoData, snp.ids=snp.ids, anom=anoms,
  centromere=centromeres.hg18)

anomStatsPlot(blData, genoData, anom.stats=stats,
  snp.ineligible=snp.failed, centromere=centromeres.hg18)

close(blData)
close(genoData)
}

\keyword{manip}
\keyword{hplot}

