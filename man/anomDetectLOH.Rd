\name{anomDetectLOH}
\alias{anomDetectLOH}
\title{
LOH Method for Chromosome Anomaly Detection 
}
\description{
\code{anomDetectLOH} breaks a chromosome up into segments of homozygous runs
of SNP markers determined by change points in Log R Ratio and 
selects segments which are likely to be anomalous.
}

\usage{
anomDetectLOH(intenData, genoData, scan.ids, chrom.ids, snp.ids,
  known.anoms, smooth = 50, min.width = 5, nperm = 10000, alpha = 0.001,
  run.size = 50, inter.size = 4, homodel.min.num = 10, homodel.thresh = 10,
  small.num = 20, small.thresh = 2.25, medium.num = 50, medium.thresh = 2,
  long.num = 100, long.thresh = 1.5, small.na.thresh = 2.5,
  length.factor = 5, merge.fac = 0.85, min.lrr.num = 20, verbose = TRUE)
}

\arguments{
  \item{intenData}{
  An \code{\link{IntensityData}} object containing the Log R Ratio.  
  The order of the rows of intenData and the snp annotation
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
  vector of (unique) chromosomes to process.  Should correspond to
  integer chromosome codes in \code{intenData}.  Recommended for use with
  autosomes, X (males will be ignored), and the
  pseudoautosomal (XY) region.
}
  \item{snp.ids}{
  vector of eligible snp ids.  Usually exclude failed and intensity-only snps.
  Also recommended to exclude an HLA region on chromosome 6 and
  XTR region on X chromosome.  See \code{\link{HLA}} and
  \code{\link{pseudoautosomal}}.
  If there are SNPs annotated in the centromere gap, exclude these as
  well (see \code{\link{centromeres}}).
}
\item{known.anoms}{
  data.frame of known anomalies (usually from \code{\link{anomDetectBAF}}); 
  must have "scanID","chromosome","left.index","right.index".
  Here "left.index" and "right.index" are row indices of intenData. Left and right
refer to start and end of anomaly, respectively, in position order.
}
\item{smooth}{
  number of markers for smoothing region.  See \code{\link{smooth.CNA}} 
  in the \pkg{\link{DNAcopy}} package.
}
\item{min.width}{
  minimum number of markers for segmenting.  See \code{\link{segment}}
  in the \pkg{\link{DNAcopy}} package.
}
\item{nperm}{
  number of permutations.
See \code{\link{segment}} in the \pkg{\link{DNAcopy}} package.
}
\item{alpha}{
  significance level. See \code{\link{segment}} in the
  \pkg{\link{DNAcopy}} package.
}

\item{run.size}{
  number of markers to declare a 'homozygous' run (here 'homozygous' 
   includes homozygous and missing) 
}
  \item{inter.size}{
number of consecutive heterozygous markers allowed to interrupt a 'homozygous' run
}
\item{homodel.min.num}{
  minimum number of markers to detect extreme difference in lrr (for homozygous deletion)
}
\item{homodel.thresh}{
  threshold for measure of deviation from non-anomalous needed to declare segment a homozygous deletion.
}
\item{small.num}{
  minimum number of SNP markers to declare segment as an anomaly 
   (other than homozygous deletion)
}
\item{small.thresh}{
 threshold for measure of deviation from non-anomalous to declare segment anomalous if
 number of SNP markers is between \code{small.num} and \code{medium.num}.
}
\item{medium.num}{
  threshold for number of SNP markers to identify 'medium' size segment 
}
\item{medium.thresh}{
  threshold for measure of deviation from non-anomalous needed to declare segment anomalous if
 number of SNP markers is between \code{medium.num} and \code{long.num}.
}
\item{long.num}{
  threshold for number of SNP markers to identify 'long' size segment
}
\item{long.thresh}{
  threshold for measure of deviation from non-anomalous when number of markers is bigger than \code{long.num}
}
\item{small.na.thresh}{
   threshold measure of deviation from non-anomalous when number of markers is between \code{small.num} and \code{medium.num}
    and 'local mad.fac' is NA.  See Details section for definition of
  'local mad.fac'.
}
\item{length.factor}{
  window around anomaly defined as \code{length.factor}*(no. of markers in segment) 
   on either side of the given segment.  Used in determining 'local mad.fac'. See Details section.
}
\item{merge.fac}{
  threshold for 'sd.fac'= number of baseline standard deviations of segment mean from baseline mean; 
  consecutive segments with 'sd.fac' above threshold are merged 
}
\item{min.lrr.num}{
  if any 'non-anomalous' interval has fewer markers than \code{min.lrr.num},
  interval is ignored in finding non-anomalous baseline unless it's the only piece left
}
\item{verbose}{
   logical indicator whether to print the scan id currently being processed
}
}

\details{
Detection of anomalies with loss of heterozygosity accompanied by change in Log R Ratio.
 Male samples for X chromosome are not processed.

Circular binary segmentation (CBS) (using the R-package \pkg{\link{DNAcopy}})
is applied to LRR values and, in parallel, runs of homozygous or missing genotypes 
of a certain minimal size (\code{run.size}) (and allowing for some interruptions 
by no more than \code{inter.size} heterozygous SNPs ) are identified.  Intervals from
\code{known.anoms} are excluded from the identification of runs.
 After some possible merging of consecutive CBS segments 
(based on satisfying a threshold  \code{merge.fac} for deviation 
from non-anomalous baseline), the homozygous runs are intersected 
with the segments from CBS. 


Determination of anomalous segments is based on 
a combination of number-of-marker thresholds and deviation from a non-anomalous 
baseline.  Segments are declared anomalous if deviation from non-anomalous is above
 corresponding thresholds. (See \code{small.num}, \code{small.thresh}, \code{medium.num},\code{medium.thresh},
\code{long.num},\code{long.thresh},and \code{small.na.thresh}.) 
 Non-anomalous median and MAD are defined for each sample-chromosome combination.
 Intervals from \code{known.anoms} and the homozygous runs 
identified are excluded; remaining regions are the non-anomalous baseline. 


Deviation from non-anomalous is measured by  
a combination of a chromosome-wide 'mad.fac' and a 'local mad.fac' (both the average
and the minimum of these two measures are used). 
 Here 'mad.fac' is (segment median-non-anomalous median)/(non-anomalous MAD) and
'local mad.fac' is the same definition except the non-anomalous median and MAD
are computed over a window including the segment (see \code{length.factor}).
Median and MADare found for eligible LRR values.
}

\value{
A list with the following elements:
\item{raw}{raw homozygous run data, not including any regions present in \code{known.anoms}.  
A data.frame with the following
  columns: Left and right
refer to start and end of anomaly, respectively, in position order.

  \itemize{
    \item\code{left.index}:  row index of intenData indicating left endpoint of segment
    \item\code{right.index}:  row index of intenData indicating right endpoint of segment
    \item\code{left.base}:  base position of left endpoint of segment
   \item\code{right.base}:  base position of right endpoint of segment
    \item \code{scanID}:  integer id of scan
    \item \code{chromosome}:  chromosome as integer code
  }
}
\item{raw.adjusted}{  data.frame of runs after merging and intersecting
  with CBS segments,  with the following
  columns: Left and right
refer to start and end of anomaly, respectively, in position order.

  \itemize{
    \item \code{scanID}:  integer id of scan
    \item \code{chromosome}:  chromosome as integer code
    \item\code{left.index}:  row index of intenData indicating left endpoint of segment
    \item\code{right.index}:  row index of intenData indicating right endpoint of segment
    \item\code{left.base}:  base position of left endpoint of segment
   \item\code{right.base}:  base position of right endpoint of segment
    \item \code{num.mark}:  number of eligible SNP markers in segment
    \item \code{seg.median}:  median of eligible LRR values in segment
    \item \code{seg.mean}:  mean of eligible LRR values in segment
    \item \code{mad.fac}:  measure of deviation from non-anomalous baseline,  equal to
         abs(median of segment - baseline median)/(baseline MAD);
         used in determining anomalous segments
    \item \code{sd.fac}:  measure of deviation from non-anomalous baseline, equal to
         abs(mean of segment - baseline mean)/(baseline standard deviation);
         used in determining whether to merge
    \item \code{local}:  measure of deviation from non-anomalous baseline used equal to
         abs(median of segment - local baseline median)/(local baseline MAD);
         local baseline consists of eligible LRR values in a window around segment;
         used in determining anomalous segments
    \item \code{num.segs}:  number of segments found by CBS for the given chromosome
    \item \code{chrom.nonanom.mad}:  MAD of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{chrom.nonanom.median}:  median of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{chrom.nonanom.mean}:  mean of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{chrom.nonanom.sd}:  standard deviation of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{sex}:  sex of the scan id coded as "M" or "F"
  }
}
\item{filtered}{ data.frame of the segments identified as anomalies.  Columns are the
  same as in \code{raw.adjusted}.
}
\item{base.info}{ data.frame with columns:
  \itemize{
    \item \code{chrom.nonanom.mad}: MAD of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{chrom.nonanom.median}: median of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{chrom.nonanom.mean}: mean of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{chrom.nonanom.sd}: standard deviation of eligible LRR values in non-anomalous regions across the chromosome
    \item \code{sex}: sex of the scan id coded as "M" or "F"
    \item \code{num.runs}: number of original homozygous runs found for given scan/chromosome
    \item \code{num.segs}: number of segments for given scan/chromosome produced by CBS
    \item \code{scanID}: integer id of scan
    \item \code{chromosome}: chromosome as integer code
    \item \code{sex}: sex of the scan id coded as "M" or "F"
  }
}
\item{segments}{ data.frame of the segmentation found by CBS with columns:
  \itemize{
    \item \code{scanID}:  integer id of scan
    \item \code{chromosome}:  chromosome as integer code
    \item\code{left.index}:  row index of intenData indicating left endpoint of segment
    \item\code{right.index}:  row index of intenData indicating right endpoint of segment
    \item\code{left.base}:  base position of left endpoint of segment
   \item\code{right.base}:  base position of right endpoint of segment
    \item \code{num.mark}:  number of eligible SNP markers in the segment
    \item \code{seg.mean}:  mean of eligible LRR values in the segment
    \item \code{sd.fac}:  measure of deviation from baseline equal to
         abs(mean of segment - baseline mean)/(baseline standard deviation)
         where the baseline is over non-anomalous regions
  }
}
\item{merge}{data.frame of scan id/chromosome pairs for which merging occurred.
  \itemize{
    \item \code{scanID}:  integer id of scan
    \item \code{chromosome}:  chromosome as integer code
  }
}
}
\references{
See references in \code{\link{segment}} in the package \pkg{\link{DNAcopy}}.
}
\author{
  Cecelia Laurie
}

\seealso{
 \code{\link{segment}} and \code{\link{smooth.CNA}} in the package \pkg{\link{DNAcopy}}, 
  also \code{\link{findBAFvariance}}, \code{\link{anomDetectLOH}}
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

# example for known.anoms, would get this from anomDetectBAF
known.anoms <- data.frame("scanID"=scan.ids[1],"chromosome"=21,
  "left.index"=100,"right.index"=200)

LOH.anom <- anomDetectLOH(blData, genoData, scan.ids=scan.ids,
  chrom.ids=chrom.ids, snp.ids=snp.ids, known.anoms=known.anoms)

close(blData)
close(genoData)
}
\keyword{manip}
