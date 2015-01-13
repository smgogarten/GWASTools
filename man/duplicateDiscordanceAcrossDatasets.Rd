\name{duplicateDiscordanceAcrossDatasets} % what should this value be? 'duplicateCalcsAcrossDatasets'..?
\alias{duplicateDiscordanceAcrossDatasets}
\alias{minorAlleleDetectionAccuracy}
\alias{dupDosageCorAcrossDatasets}
\title{Functions to check discordance and allelic dosage correlation across datasets}

\description{
  These functions compare genotypes in pairs of duplicate scans of the same sample
  across multiple datasets. 'duplicateDiscordanceAcrossDatasets' finds the number of discordant genotypes both by scan and by SNP.
  'dupDosageCorAcrossDatasets' calculates squared correlation (r2) between allelic dosages both by scan and by SNP, allowing for
  comparision between imputed datasets or between imputed and observed -
  i.e., where one or more of the datasets contains continuous dosage [0,2] rather than discrete allele counts \{0,1,2\}.
}

\usage{
duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
  match.snps.on=c("position", "alleles"),
  subjName.cols, snpName.cols=NULL,
  one.pair.per.subj=TRUE, minor.allele.only=FALSE,
  missing.fail=c(FALSE, FALSE),
  scan.exclude1=NULL, scan.exclude2=NULL,
  snp.exclude1=NULL, snp.exclude2=NULL,
  snp.include=NULL,
  verbose=TRUE)

minorAlleleDetectionAccuracy(genoData1, genoData2,
  match.snps.on=c("position", "alleles"),
  subjName.cols, snpName.cols=NULL,
  missing.fail=TRUE,
  scan.exclude1=NULL, scan.exclude2=NULL,
  snp.exclude1=NULL, snp.exclude2=NULL,
  snp.include=NULL,
  verbose=TRUE)

dupDosageCorAcrossDatasets(genoData1, genoData2,
  match.snps.on=c("position", "alleles"),
  subjName.cols="subjectID", snpName.cols=NULL,
  scan.exclude1=NULL, scan.exclude2=NULL,
  snp.exclude1=NULL, snp.exclude2=NULL,
  snp.include=NULL,
  snp.block.size=5000, scan.block.size=100,
  verbose=TRUE)
}

\arguments{
  \item{genoData1}{\code{\link{GenotypeData}} object containing the
    first dataset.}
  \item{genoData2}{\code{\link{GenotypeData}} object containing the
    second dataset.}
  \item{match.snps.on}{One or more of ("position", "alleles", "name") indicating how to match SNPs.  "position" will match SNPs on chromosome and position, "alleles" will also require the same alleles (but A/B designations need not be the same), and "name" will match on the columns give in \code{snpName.cols}.}
  \item{subjName.cols}{2-element character vector indicating the names of the
    annotation variables that will be identical for duplicate scans in
    the two datasets.  (Alternatively, one character value that will be recycled).}
  \item{snpName.cols}{2-element character vector indicating the names of the
    annotation variables that will be identical for the same SNPs in the
    two datasets.  (Alternatively, one character value that will be recycled).}
\item{one.pair.per.subj}{A logical indicating whether a single pair of
  scans should be randomly selected for each subject with more than 2 scans.}
\item{minor.allele.only}{A logical indicating whether discordance should be calculated
  only between pairs of scans in which at least one scan has a genotype
  with the minor allele (i.e., exclude major allele homozygotes).}
\item{missing.fail}{For \code{duplicateDiscordanceAcrossDatasets},
    a 2-element logical vector indicating whether missing
  values in datasets 1 and 2, respectively, will be considered failures
  (discordances with called genotypes in the other dataset).  For
  \code{minorAlleleDetectionAccuracy}, a single logical indicating
  whether missing values in dataset 2 will be considered false negatives
  (\code{missing.fail=TRUE}) or ignored (\code{missing.fail=FALSE}).}
\item{scan.exclude1}{An integer vector containing the ids of scans to be
  excluded from the first dataset. }
\item{scan.exclude2}{An integer vector containing the ids of scans to be
  excluded from the second dataset. }
\item{snp.exclude1}{An integer vector containing the ids of SNPs to be
  excluded from the first dataset. }
\item{snp.exclude2}{An integer vector containing the ids of SNPs to be
  excluded from the second dataset. }
  \item{snp.include}{List of SNPs to include in the comparison.  Should
    match the contents of the columns referred to by \code{snpName.cols}.
    Only valid if \code{match.snps.on} includes "name".}
\item{snp.block.size}{Block size for SNPs}
\item{scan.block.size}{Block size for scans}
\item{verbose}{Logical value specifying whether to show progress
  information.}
}

\details{
  \code{duplicateDiscordanceAcrossDatasets} calculates discordance metrics both by
  scan and by SNP.  If \code{one.pair.per.subj=TRUE} (the default), each
  subject with more than two duplicate genotyping instances will have
  one scan from each dataset randomly selected for computing discordance.  If
  \code{one.pair.per.subj=FALSE}, discordances will be calculated
  pair-wise for all possible cross-dataset pairs for each subject.

  \code{dupDosageCorAcrossDatasets} calculates squared dosage correlation both
by scan and by SNP. Note it only allows for one pair of duplicate scans
per sample.  For this function only, \code{genoData1} and
\code{genoData2} must have been created with
\code{\link{GdsGenotypeReader}} objects.

By default, overlapping variants are identified based on position and
alleles. Alleles are determined via 'getAlleleA' and 'getAlleleB'
accessors, so users should ensure these variables are referring to the
same strand orientation in both datests (e.g., both plus strand
alleles). It is not necessary for the A/B ordering to be consistent
across datasets. For example, two variants at the same position with
alleleA="C" and alleleB="T" in genoData1 and alleleA="T" and alleleB="C"
in genoData2 will stil be identified as overlapping.

  If \code{minor.allele.only=TRUE}, the allele frequency will be
  calculated in \code{genoData1}, using only samples common to both datasets.

  If \code{snp.include=NULL} (the default), discordances will be found
  for all SNPs common to both datasets.

  \code{genoData1} and \code{genoData2} should each have "alleleA" and
  "alleleB" defined in their SNP annotation.  If allele coding cannot be
  found, the two datasets are assumed to have identical coding. Note
  that 'dupDosageCorAcrossDatasets' can NOT detect where strand-ambiguous (A/T or
C/G) SNPs are annotated on different strands, although the r2 in these
instances would be unaffected: r may be negative but r2 will be positive.

  \code{minorAlleleDetectionAccuracy} summarizes the accuracy of minor
  allele detection in \code{genoData2} with respect to \code{genoData1}
  (the "gold standard").
  \code{TP=}number of true positives, \code{TN=}number of true negatives,
  \code{FP=}number of false positives, and \code{FN=}number of false
  negatives.
  Accuracy is represented by four metrics:
  \itemize{
    \item sensitivity for each SNP as \code{TP/(TP+FN)}
    \item specificity for each SNP as \code{TN/(TN+FP)}
    \item positive predictive value for each SNP as \code{TP/(TP+FP)}
    \item negative predictive value for each SNP as \code{TN/(TN+FN)}.
  }
  \code{TP}, \code{TN}, \code{FP}, and \code{FN} are calculated as follows:
  \tabular{lllll}{
              \tab    \tab           \tab genoData1 \tab \cr
              \tab    \tab mm        \tab Mm        \tab MM\cr
              \tab mm \tab 2TP       \tab 1TP + 1FP \tab 2FP\cr
    genoData2 \tab Mm \tab 1TP + 1FN \tab 1TN + 1TP \tab 1TN + 1FP\cr
	      \tab MM \tab 2FN       \tab 1FN + 1TN \tab 2TN\cr
              \tab -- \tab 2FN       \tab 1FN       \tab \cr
	    }
  "M" is the major allele and
  "m" is the minor allele (as calculated in \code{genoData1}).
  "-" is a missing call in \code{genoData2}.
  Missing calls in \code{genoData1} are ignored. If
  \code{missing.fail=FALSE}, missing calls in \code{genoData2} (the last
  row of the table) are also ignored.
}

\value{
  SNP annotation columns returned by all functions are:
  \item{chromosome}{chromosome}
  \item{position}{base pair position}
  \item{snpID1}{snpID from genoData1}
  \item{snpID2}{snpID from genoData2}

  If matching on "alleles":
  \item{alleles}{alleles sorted alphabetically}
  \item{alleleA1}{allele A from genoData1}
  \item{alleleB1}{allele B from genoData2}
  \item{alleleA2}{allele A from genoData2}
  \item{alleleB2}{allele B from genoData2}

  If matching on "name":
  \item{name}{the common SNP name given in \code{snpName.cols}}

  \code{duplicateDiscordanceAcrossDatasets} returns a list with two data
  frames:
  The data.frame "discordance.by.snp" contains the
  SNP annotation columns listed above as well as:
  \item{discordant}{number of discordant pairs}
  \item{npair}{number of pairs examined}
  \item{n.disc.subj}{number of subjects with at least one discordance}
  \item{discord.rate}{discordance rate i.e. discordant/npair}

  The data.frame "discordance.by.subject" contains a list of matrices (one for each subject) with the pair-wise discordance between the different genotyping instances of the subject.


  \code{minorAlleleDetectionAccuracy} returns a data.frame with the
  SNP annotation columns listed above as well as:
  \item{npair}{number of sample pairs compared (non-missing in \code{genoData1})}
  \item{sensitivity}{sensitivity}
  \item{specificity}{specificity}
  \item{positivePredictiveValue}{Positive predictive value}
  \item{negativePredictiveValue}{Negative predictive value}


  \code{dupDosageCorAcrossDatasets} returns a list with two data
  frames:

  The data.frame "snps" contains the by-SNP r2 values with the
  SNP annotation columns listed above as well as:
  \item{nsamp.dosageR2}{number of samples in r2 calculation (i.e., non
  missing data in both genoData1 and genoData2)}
  \item{dosageR2}{squared dosage correlation}

    The data.frame "samps" contains the by-sample r2 values with the
    following columns:
    \item{subjectID}{subject-level identifier for duplicate sample pair}
    \item{scanID1}{scanID from genoData1}
    \item{scanID2}{scanID from genoData2}
  \item{nsnp.dosageR2}{number of SNPs in r2 calculation (i.e., non
  missing data in both genoData1 and genoData2)}
    \item{dosageR2}{squared dosage correlation}

  If no duplicate scans or no common SNPs are found, these functions issue a warning
  message and return \code{NULL}.
}

\author{Stephanie Gogarten, Jess Shen, Sarah Nelson}

\seealso{\code{\link{GenotypeData}}, \code{\link{duplicateDiscordance}},
  \code{\link{duplicateDiscordanceProbability}}}

\examples{
# first set
snp1 <- data.frame(snpID=1:10, chromosome=1L, position=101:110,
                   rsID=paste("rs", 101:110, sep=""),
                   alleleA="A", alleleB="G", stringsAsFactors=FALSE)
scan1 <- data.frame(scanID=1:3, subjectID=c("A","B","C"), sex="F", stringsAsFactors=FALSE)
mgr <- MatrixGenotypeReader(genotype=matrix(c(0,1,2), ncol=3, nrow=10), snpID=snp1$snpID,
                            chromosome=snp1$chromosome, position=snp1$position, scanID=1:3)
genoData1 <- GenotypeData(mgr, snpAnnot=SnpAnnotationDataFrame(snp1),
                          scanAnnot=ScanAnnotationDataFrame(scan1))

# second set
snp2 <- data.frame(snpID=1:5, chromosome=1L,
                   position=as.integer(c(101,103,105,107,107)),
                   rsID=c("rs101", "rs103", "rs105", "rs107", "rsXXX"),
                   alleleA= c("A","C","G","A","A"),
                   alleleB=c("G","T","A","G","G"),
                   stringsAsFactors=FALSE)
scan2 <- data.frame(scanID=1:3, subjectID=c("A","C","C"), sex="F", stringsAsFactors=FALSE)
mgr <- MatrixGenotypeReader(genotype=matrix(c(1,2,0), ncol=3, nrow=5), snpID=snp2$snpID,
                            chromosome=snp2$chromosome, position=snp2$position, scanID=1:3)
genoData2 <- GenotypeData(mgr, snpAnnot=SnpAnnotationDataFrame(snp2),
                          scanAnnot=ScanAnnotationDataFrame(scan2))

duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
  match.snps.on="position",
  subjName.cols="subjectID")

duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
  match.snps.on=c("position", "alleles"),
  subjName.cols="subjectID")

duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
  match.snps.on=c("position", "alleles", "name"),
  subjName.cols="subjectID",
  snpName.cols="rsID")

duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
  subjName.cols="subjectID",
  one.pair.per.subj=FALSE)

minorAlleleDetectionAccuracy(genoData1, genoData2,
  subjName.cols="subjectID")

dupDosageCorAcrossDatasets(genoData1, genoData2,
  scan.exclude2=scan2$scanID[duplicated(scan2$subjectID)])
}

\keyword{manip}
