\name{anomDetectBAF}
\alias{anomSegmentBAF}
\alias{anomFilterBAF}
\alias{anomDetectBAF}
\title{
BAF Method for Chromosome Anomaly Detection 
}
\description{
\code{anomSegmentBAF} for each sample and chromosome, breaks the chromosome up into 
segments marked by change points of a metric based on B Allele Frequency (BAF) values. 


\code{anomFilterBAF} selects segments which are likely to be anomalous.

\code{anomDetectBAF} is a wrapper to run \code{anomSegmentBAF} and
\code{anomFilterBAF} in one step.
}

\usage{
anomSegmentBAF(intenData, genoData, scan.ids, chrom.ids, snp.ids,
  smooth = 50, min.width = 5, nperm = 10000, alpha = 0.001,
  verbose = TRUE)

anomFilterBAF(intenData, genoData, segments, snp.ids, centromere,
  low.qual.ids = NULL, num.mark.thresh = 15, long.num.mark.thresh = 200,
  sd.reg = 2, sd.long = 1, low.frac.used = 0.1, run.size = 10,
  inter.size = 2, low.frac.used.num.mark = 30, very.low.frac.used = 0.01, 
  low.qual.frac.num.mark = 150, lrr.cut = -2, ct.thresh = 10,
  frac.thresh = 0.1, verbose=TRUE)

anomDetectBAF(intenData, genoData, scan.ids, chrom.ids, snp.ids,
  centromere, low.qual.ids = NULL, ...)
}
\arguments{
\item{intenData}{
  An \code{\link{IntensityData}} object containing the B Allele
  Frequency.  The order of the rows of intenData and the snp annotation
  are expected to be by chromosome and then by position within chromosome.
  The scan annotation should contain sex, coded as "M" for male and
  "F" for female.
}
\item{genoData}{
  A \code{\link{GenotypeData}} object.  The order of the rows of genoData 
  and the snp annotation are expected to be by chromosome and then 
  by position within chromosome.
}
\item{scan.ids}{
  vector of scan ids (sample numbers) to process
}
\item{chrom.ids}{
  vector of (unique) chromosomes to process.  Recommended to include
  all autosomes.
}
\item{snp.ids}{
  vector of eligible snp ids.  Usually exclude failed and intensity-only SNPs.
  Also recommended to exclude an HLA region on chromosome 6 and
  XTR region on chromosome 23 (X).  See \code{\link{HLA}} and \code{\link{pseudoautosomal}}.
}
\item{smooth}{
  number of markers for smoothing region.  See \code{\link{smooth.CNA}} 
  in the \code{\link{DNAcopy}} package.
}
\item{min.width}{
  minimum number of markers for a segment.  See \code{\link{segment}}
  in the \code{\link{DNAcopy}} package.
}
\item{nperm}{
  number of permutations for deciding significance in segmentation.
   See \code{\link{segment}} in the \code{\link{DNAcopy}} package.
}
\item{alpha}{
  significance level.  See \code{\link{segment}} in the
  \code{\link{DNAcopy}} package.
}
\item{verbose}{
   logical indicator whether to print information about the scan id currently being processed.
 anomSegmentBAF prints each scan id; anomFilterBAF prints a message after every 10 samples: "processing ith scan id out of n"
  where "ith" with be 10, 10, etc. and "n" is the total number of samples
}
\item{segments}{
  data.frame of segments from \code{anomSegmentBAF}.  Names must
  include "scanID", "chromosome", "num.mark", "left.index", "right.index", "seg.mean".
 Here "left.index" and "right.index" are row indices of intenData. Left and right
refer to start and end of anomaly,respectively, in position order.
}
\item{centromere}{
  data.frame with centromere position information. Names must include
  "chrom", "left.base", "right.base".  Valid values for "chrom" are
  1:22, "X", "Y", "XY".  Here "left.base" and "right.base"
  are base positions of start and end of centromere location in position order.

}
\item{low.qual.ids}{
  scan ids determined to be low quality for which some segments are filtered
  based on more stringent criteria.  Default is NULL.  Usual choice are
  scan ids for which median BAF across autosomes > 0.05.  See
  \code{\link{sdByScanChromWindow}} and \code{\link{medianSdOverAutosomes}}.
}
\item{num.mark.thresh}{
  minimum number of SNP markers in a segment to be considered for anomaly
}
\item{long.num.mark.thresh}{
  min number of markers for "long" segment to be considered for anomaly
  for which significance threshold criterion is allowed to be less stringent
}
\item{sd.reg}{
  number of baseline standard deviations of segment mean from a baseline
  mean for "normal" needed to declare segment anomalous. This number is given by 
 abs(mean of segment - baseline mean)/(baseline standard deviation)
}
\item{sd.long}{
  same meaning as \code{sd.reg} but applied to "long" segments
}
\item{low.frac.used}{
  if fraction of heterozygous or missing SNP markers compared with number of 
  eligible SNP markers in segment is below this, more stringent criteria 
  are applied to declare them anomalous. 
}
\item{run.size}{
  min length of run of missing or heterozygous SNP markers for possible 
  determination of homozygous deletions 
}
\item{inter.size}{
  number of homozygotes allowed to "interrupt" run for possible
  determination of homozygous deletions 
}
\item{low.frac.used.num.mark}{
  number of markers threshold for \code{low.frac.used} segments (which are not
  declared homozygous deletions 
}
\item{very.low.frac.used}{
  any segments with (num.mark)/(number of markers in interval) less than this 
  are filtered out since they tend to be false positives
}
\item{low.qual.frac.num.mark}{
  minimum num.mark threshold for low quality scans (\code{low.qual.ids}) 
  for segments that are also below low.frac.used threshold
}
\item{lrr.cut}{
  look for runs of LRR values below \code{lrr.cut} to adjust homozygous deletion endpoints
}
\item{ct.thresh}{
  minimum number of LRR values below \code{lrr.cut} needed in order to adjust
}
\item{frac.thresh}{
  investigate interval for homozygous deletion only if \code{lrr.cut} and \code{ct.thresh}
  thresholds met and (# LRR values below \code{lrr.cut})/(# eligible SNPs in segment) > \code{frac.thresh}
}
\item{...}{
  arguments to pass to \code{anomFilterBAF}
}
}
\details{
\code{anomSegmentBAF} uses the function \code{\link{segment}} from
the \code{DNAcopy} package to perform circular binary segmentation
on a metric based on BAF values.  The metric for a given sample/chromosome 
is sqrt(min(BAF,1-BAF,abs(BAF-median(BAF))) where the median is 
across BAF values on the chromosome.  Only BAF values for heterozygous or
missing SNPs are used.

\code{anomFilterBAF} determines anomalous segments based on a combination
of thresholds for number of SNP markers in the segment and on deviation from
a "normal" baseline.  (See \code{num.mark.thresh},\code{long.num.mark.thresh},
\code{sd.reg}, and \code{sd.long}.)  The "normal" baseline metric mean and standard deviation
are found across all autosomes not segmented by \code{anomSegmentBAF}.  This is why
it is recommended to include all autosomes for the argument \code{chrom.ids} to
ensure a more accurate baseline. 
 

Some initial filtering is done,
including possible merging of consecutive segments meeting \code{sd.reg}
threshold along with other criteria (such as not spanning the centromere)
 and adjustment for accurate
break points for possible homozygous deletions (see \code{lrr.cut},
 \code{ct.thresh}, \code{frac.thresh}, \code{run.size}, and \code{inter.size}).
 Male samples for chromosome 23 (X) are not processed.

 More stringent criteria are applied to some segments 
(see \code{low.frac.used},\code{low.frac.used.num.mark}, 
\code{very.low.frac.used}, \code{low.qual.ids}, and
\code{low.qual.frac.num.mark}).

\code{anomDetectBAF} runs \code{anomSegmentBAF} with default values and
then runs \code{anomFilterBAF}.  Additional parameters for
\code{anomFilterBAF} may be passed as arguments.
}

\value{
\code{anomSegmentBAF} returns a data.frame with the following elements: Left and right
refer to start and end of anomaly, respectively, in position order.

\item{scanID}{integer id of scan}
\item{chromosome}{chromosome as integer where 23 refers to X chromosome}
\item{left.index}{row index of intenData indicating left endpoint of segment}
\item{right.index}{row index of intenData indicating right endpoint of segment}
\item{num.mark}{number of heterozygous or missing SNPs in the segment }
\item{seg.mean}{mean of the BAF metric over the segment

}

   


\code{anomFilterBAF} and \code{anomDetectBAF} return a list with the
following elements:
\item{raw}{data.frame of raw segmentation data, with same output as
  \code{anomSegmentBAF} as well as:
   \itemize{
     \item \code{left.base}:  base position of left endpoint of segment
      \item \code{right.base}:  base position of right endpoint of segment
      \item \code{sex}:  sex of scan.id coded as "M" or "F"
      \item \code{sd.fac}:  measure of deviation from baseline equal to
         abs(mean of segment - baseline mean)/(baseline standard deviation);
         used in determining anomalous segments
   }
}   
\item{filtered}{data.frame of the segments identified as anomalies, with the same columns as
  \code{raw} as well as:
  \itemize{
    \item \code{merge}:  TRUE if segment was a result of merging. Consecutive segments 
       from output of \code{anomSegmentBAF} that meet certain criteria are merged. 
    \item \code{homodel.adjust}:  TRUE if original segment was adjusted to
    narrow in on a homozygous deletion
    \item \code{frac.used}:  fraction of (eligible) heterozygous or missing SNP markers compared with total number of 
  eligible SNP markers in segment
  }
}
\item{base.info}{data frame with columns:
  \itemize{
    \item \code{scanID}:  integer id of scan
    \item \code{base.mean}:  mean of non-anomalous baseline.  This is the mean of the
      BAF metric for heterozygous and missing SNPs over all unsegmented autosomes
      that were considered.
    \item \code{base.sd}:  standard deviation of non-anomalous baseline
    \item \code{chr.ct}:  number of unsegmented chromosomes used in determining
    the non-anomalous baseline
  }
}
\item{seg.info}{data frame with columns:
  \itemize{
    \item \code{scanID}: integer id of scan
    \item \code{chromosome}: chromosome as integer
    \item \code{num.segs}: number of segments produced by \code{anomSegmentBAF}
  }
}
}
\references{
%% ~put references to the literature/web site here ~
See references in \code{\link{segment}} in the package \code{\link{DNAcopy}}.
 The BAF metric used is modified from Itsara,A., \emph{et.al} (2009) Population
Analysis of Large Copy Number Variants and Hotspots of Human Genetic Disease.
\emph{American Journal of Human Genetics}, \bold{84}, 148--161.
}
\author{
Cecelia Laurie
}

\note{
It is recommended to include all autosomes as input.  This ensures a more
accurate determination of baseline information.
}

\seealso{
  \code{\link{segment}} and \code{\link{smooth.CNA}} in the package \code{\link{DNAcopy}}, 
  also \code{\link{findBAFvariance}}, \code{\link{anomDetectLOH}}
}
\examples{
library(GWASdata)
data(illuminaScanADF)
data(illuminaSnpADF)

blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
blnc <- NcdfIntensityReader(blfile)
blData <-  IntensityData(blnc, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

genofile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
genonc <- NcdfGenotypeReader(genofile)
genoData <-  GenotypeData(genonc, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

# segment BAF
scan.ids <- illuminaScanADF$scanID[1:2]
chrom.ids <- unique(illuminaSnpADF$chromosome)
snp.ids <- illuminaSnpADF$snpID[illuminaSnpADF$missing.n1 < 1]
seg <- anomSegmentBAF(blData, genoData, scan.ids=scan.ids,
                      chrom.ids=chrom.ids, snp.ids=snp.ids)

# filter segments to detect anomalies
data(centromeres.hg18)
filt <- anomFilterBAF(blData, genoData, segments=seg, snp.ids=snp.ids,
                      centromere=centromeres.hg18)

# alternatively, run both steps at once
anom <- anomDetectBAF(blData, genoData, scan.ids=scan.ids, chrom.ids=chrom.ids,
                      snp.ids=snp.ids, centromere=centromeres.hg18)
}
\keyword{manip}
