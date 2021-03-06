\name{batchTest}
\alias{batchChisqTest}
\alias{batchFisherTest}
\title{Batch Effects of Genotyping}
\description{
  \code{batchChisqTest} calculates Chi-square values for batches from
  2-by-2 tables of SNPs, comparing each batch with the other batches.
  \code{batchFisherTest} calculates Fisher's exact test values.
}

\usage{
batchChisqTest(genoData, batchVar, snp.include = NULL,
               chrom.include = 1:22, sex.include = c("M", "F"),
               scan.exclude = NULL, return.by.snp = FALSE,
               correct = TRUE, verbose = TRUE)

batchFisherTest(genoData, batchVar, snp.include = NULL,
                chrom.include = 1:22, sex.include = c("M", "F"),
                scan.exclude = NULL, return.by.snp = FALSE,
                conf.int = FALSE, verbose = TRUE)
}
		
\arguments{
  \item{genoData}{\code{\link{GenotypeData}} object}
  \item{batchVar}{A character string indicating which annotation
    variable should be used as the batch.}
  \item{snp.include}{A vector containing the IDs of SNPs to include.} 
  \item{chrom.include}{Integer vector with codes for chromosomes to
    include. Ignored if \code{snp.include} is not \code{NULL}. Default is 1:22 (autosomes).  Use 23, 24, 25, 26, 27 for X, XY, Y, M, Unmapped
    respectively}
  \item{sex.include}{Character vector with sex to include.  Default is
    c("M", "F").  If sex chromosomes are present in \code{chrom.include},
    only one sex is allowed.}
  \item{scan.exclude}{A vector containing the IDs of scans to
    be excluded.}
  \item{return.by.snp}{Logical value to indicate whether snp-by-batch
  matrices should be returned.}
  \item{conf.int}{Logical value to indicate if a confidence interval should be computed.}
  \item{correct}{Logical value to specify whether to apply the Yates continuity correction.}
  \item{verbose}{Logical value specifying whether to show progress information.}
}

\details{	
  Because of potential batch effects due to sample processing and genotype calling, batches are an important experimental design factor.
  
  \code{batchChisqTest} calculates the Chi square values from 2-by-2
  table for each SNP, comparing each batch with the other batches.
  
  \code{batchFisherTest} calculates Fisher's Exact Test from 2-by-2
  table for each SNP, comparing each batch with the other batches.
  
  For each SNP and each batch, batch effect is evaluated by a 2-by-2 table:
  # of A alleles, and # of B alleles in the batch, versus
  # of A alleles, and # of B alleles in the other batches.
  Monomorphic SNPs are set to \code{NA} for all batches.
  
  The default behavior is to combine allele frequencies from males and
  females and return results for autosomes only.  If results for sex
  chromosomes (X or Y) are desired, use \code{chrom.include} with values
  23 and/or 25 and \code{sex.include}="M" or "F".

  If there are only two batches, the calculation is only performed once
  and the values for each batch will be identical.
}

\value{
  \code{batchChisqTest} returns a list with the following elements:
  \item{mean.chisq}{a vector of mean chi-squared values for each batch.}
  \item{lambda}{a vector of genomic inflation factor computed as \code{median(chisq)
      / 0.456} for each batch.}
  \item{chisq}{a matrix of chi-squared values with SNPs as rows and
    batches as columns.  Only returned if \code{return.by.snp=TRUE}.}
  
  \code{batchFisherTest} returns a list with the following elements:
  \item{mean.or}{a vector of mean odds-ratio values for each
  batch. \code{mean.or} is computed as \code{1/mean(pmin(or, 1/or))} since the odds ratio is >1 when the batch has a higher allele frequency than the other batches and <1 for the reverse.}
  \item{lambda}{a vector of genomic inflation factor computed as
    \code{median(-2*log(pval) / 1.39} for each batch.}
  Each of the following is a matrix with SNPs as rows and batches as
  columns, and is only returned if \code{return.by.snp=TRUE}:
  \item{pval}{P value}
  \item{oddsratio}{Odds ratio}
  \item{confint.low}{Low value of the confidence interval for the odds
  ratio.  Only returned if \code{conf.int=TRUE}.}
  \item{confint.high}{High value of the confidence interval for the odds
  ratio.  Only returned if \code{conf.int=TRUE}.}

  \code{batchChisqTest} and \code{batchFisherTest} both also return the
  following if \code{return.by.snp=TRUE}:
  \item{allele.counts}{matrix with total number of A and B alleles over
  all batches.}
  \item{min.exp.freq}{matrix of minimum expected allele frequency with SNPs as rows and
    batches as columns.}
}

\author{Xiuwen Zheng, Stephanie Gogarten}

\seealso{\code{\link{GenotypeData}}, \code{\link{chisq.test}},
  \code{\link{fisher.test}}
}

\examples{
library(GWASdata)
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)
data(illuminaScanADF)
genoData <-  GenotypeData(gds, scanAnnot=illuminaScanADF)

# autosomes only, sexes combined (default)
res.chisq <- batchChisqTest(genoData, batchVar="plate")
res.chisq$mean.chisq
res.chisq$lambda

# X chromosome for females
res.chisq <- batchChisqTest(genoData, batchVar="status",
  chrom.include=23, sex.include="F", return.by.snp=TRUE)
head(res.chisq$chisq)

# Fisher exact test of "status" on X chromosome for females
res.fisher <- batchFisherTest(genoData, batchVar="status",
  chrom.include=23, sex.include="F", return.by.snp=TRUE)
qqPlot(res.fisher$pval)

close(genoData)
}

\keyword{htest}
