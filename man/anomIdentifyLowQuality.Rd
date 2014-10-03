\name{anomIdentifyLowQuality}
\alias{anomIdentifyLowQuality}
\title{Identify low quality samples }

\description{
  Identify low quality samples for which false positive rate for anomaly detection 
   is likely to be high. Measures of noise (high variance) and high segmentation
  are used.
}

\usage{
anomIdentifyLowQuality(snp.annot, med.sd, seg.info,
  sd.thresh, sng.seg.thresh, auto.seg.thresh)
}

\arguments{
  \item{snp.annot}{
    \code{\link{SnpAnnotationDataFrame}} with column "eligible", where
    "eligible" is a logical vector indicating whether a SNP is eligible
    for consideration in anomaly detection (usually FALSE for HLA and
    XTR regions, failed SNPs, and intensity-only SNPs). 
   See \code{\link{HLA}} and \code{\link{pseudoautosomal}}.
  }
  \item{med.sd}{
    data.frame of median standard deviation of BAlleleFrequency (BAF) or LogRRatio (LRR) values across autosomes for each scan, with
    columns "scanID" and "med.sd".  Usually the result of
    \code{\link{medianSdOverAutosomes}}. Usually only eligible SNPs are used in
these computations. In addition, for BAF, homozygous SNPS are excluded.
  }
  \item{seg.info}{
    data.frame with segmentation information from \code{\link{anomDetectBAF}} or
    \code{\link{anomDetectLOH}}.  Columns must include "scanID",
    "chromosome", and "num.segs".  (For \code{\link{anomDetectBAF}}, segmentation information 
    is found in $seg.info from output.  For \code{\link{anomDetectLOH}}, segmentation information
     is found in $base.info from output.)
  }
  \item{sd.thresh}{ Threshold for \code{med.sd} above which scan is identified
   as low quality.
    Suggested values are 0.1 for BAF and 0.25 for LOH.
  }
  \item{sng.seg.thresh}{ Threshold for segmentation factor for a given chromosome,
  above which the chromosome is said to be highly segmented.  See Details.
    Suggested values are 0.0008 for BAF and 0.0048 for LOH.
  }
  \item{auto.seg.thresh}{Threshold for segmentation factor across autosome,
  above which the scan is said to be highly segmented.  See Details.
   Suggested values are 0.0001 for BAF and 0.0006 for LOH.
  }
}

\details{Low quality samples are determined separately with regard to each
 of the two methods of segmentation, \code{\link{anomDetectBAF}} and
 \code{\link{anomDetectLOH}}.  BAF anomalies (respectively 
 LOH anomalies) found for samples identified as low quality for BAF (respectively
 LOH) tend to have a high false positive rate.


 A scan is identified as low quality due to high variance (noise), i.e.
 if \code{med.sd} is above a certain threshold \code{sd.thresh}.


 High segmentation is often an indication of artifactual patterns in the 
  B Allele Frequency (BAF) or Log R Ratio values (LRR) that are not always captured
  by high variance.  Here segmentation information is determined by
  \code{\link{anomDetectBAF}} or  \code{\link{anomDetectLOH}} which use
 circular binary segmentation implemented by the R-package \pkg{\link{DNAcopy}}.
  The measure for high segmentation is a "segmentation factor" =
   (number of segments)/(number of eligible SNPS).  A single chromosome segmentation
factor uses information for one chromosome.  A segmentation factor across 
autosomes uses the total number of segments and eligible SNPs across all autosomes.
 See \code{med.sd}, \code{sd.thresh}, \code{sng.seg.thresh}, and \code{auto.seg.thresh}.
}

\value{
  A data.frame with the following columns:
  \item{scanID}{integer id for the scan}
  \item{chrX.num.segs}{number of segments for chromosome X}
  \item{chrX.fac}{segmentation factor for chromosome X}
  \item{max.autosome}{autosome with highest single segmentation factor}
  \item{max.auto.fac}{segmentation factor for chromosome = \code{max.autosome}}
  \item{max.auto.num.segs}{number of segments for chromosome = \code{max.autosome} }
  \item{num.ch.segd}{number of chromosomes segmented, i.e. for which change points were found}
  \item{fac.all.auto}{segmentation factor across all autosomes}
  \item{med.sd}{median standard deviation of BAF (or LRR values) across autosomes. See \code{med.sd} in Arguments section.}
  \item{type}{one of the following, indicating reason for identification as low quality:
    \itemize{
     \item\code{auto.seg}:  segmentation factor \code{fac.all.auto} above \code{auto.seg.thresh} but \code{med.sd} acceptable 
     \item \code{sd}:  standard deviation factor \code{med.sd} above \code{sd.thresh} but \code{fac.all.auto} acceptable
     \item \code{both.sd.seg}:  both high variance and high segmentation factors, \code{fac.all.auto} and \code{med.sd}, are above respective thresholds 
     \item \code{sng.seg}:  segmentation factor \code{max.auto.fac} is above \code{sng.seg.thresh} but other measures acceptable
     \item \code{sng.seg.X}: segmentation factor \code{chrX.fac} is above \code{sng.seg.thresh} but other measures acceptable
   }
}
}

\author{Cecelia Laurie}

\seealso{\code{\link{findBAFvariance}},
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

# initial scan for low quality with median SD
baf.sd <- sdByScanChromWindow(blData, genoData)
med.baf.sd <- medianSdOverAutosomes(baf.sd)
low.qual.ids <- med.baf.sd$scanID[med.baf.sd$med.sd > 0.05]

# segment and filter BAF
scan.ids <- illuminaScanADF$scanID[1:2]
chrom.ids <- unique(illuminaSnpADF$chromosome)
snp.ids <- illuminaSnpADF$snpID[illuminaSnpADF$missing.n1 < 1]
data(centromeres.hg18)
anom <- anomDetectBAF(blData, genoData, scan.ids=scan.ids, chrom.ids=chrom.ids,
  snp.ids=snp.ids, centromere=centromeres.hg18, low.qual.ids=low.qual.ids)

# further screen for low quality scans
snp.annot <- illuminaSnpADF
snp.annot$eligible <- snp.annot$missing.n1 < 1
low.qual <- anomIdentifyLowQuality(snp.annot, med.baf.sd, anom$seg.info,
  sd.thresh=0.1, sng.seg.thresh=0.0008, auto.seg.thresh=0.0001)

close(blData)
close(genoData)
}

\keyword{manip}
