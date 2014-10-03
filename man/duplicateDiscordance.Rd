\name{duplicateDiscordance}
\alias{duplicateDiscordance}
\title{Duplicate discordance}
\description{
A function to compute pair-wise genotype discordances between multiple genotyping instances of the same subject.
}
\usage{
duplicateDiscordance(genoData, subjName.col,
                     one.pair.per.subj=TRUE, corr.by.snp=FALSE,
                     minor.allele.only=FALSE, allele.freq=NULL,
                     scan.exclude=NULL, snp.exclude=NULL,
                     verbose=TRUE)
}
\arguments{
\item{genoData}{\code{\link{GenotypeData}} object} 
\item{subjName.col}{A character string indicating the name of the annotation variable
  that will be identical for duplicate scans.}
\item{one.pair.per.subj}{A logical indicating whether a single pair of
  scans should be randomly selected for each subject with more than 2 scans.}
\item{corr.by.snp}{A logical indicating whether correlation by SNP
  should be computed (may significantly increase run time).}
\item{minor.allele.only}{A logical indicating whether discordance should be calculated
  only between pairs of scans in which at least one scan has a genotype
  with the minor allele (i.e., exclude major allele homozygotes).}
\item{allele.freq}{A numeric vector with the frequency of the A allele
  for each SNP in \code{genoData}.  Required if \code{minor.allele.only=TRUE}.}
\item{scan.exclude}{An integer vector containing the ids of scans to be
  excluded. }
\item{snp.exclude}{An integer vector containing the ids of SNPs to be excluded.}
\item{verbose}{Logical value specifying whether to show progress information.} 
}
\details{
  \code{duplicateDiscordance} calculates discordance metrics both by
  scan and by SNP.  If \code{one.pair.per.subj=TRUE} (the default), each
  subject with more than two duplicate genotyping instances will have
  two scans randomly selected for computing discordance.  If
  \code{one.pair.per.subj=FALSE}, discordances will be calculated
  pair-wise for all possible pairs for each subject.
}
\value{
  A list with the following components:
  \item{discordance.by.snp}{data frame with 5  columns: 1. snpID, 2. discordant (number of discordant pairs), 3. npair (number of pairs examined), 4. n.disc.subj (number of subjects with at least one discordance), 5. discord.rate (discordance rate i.e. discordant/npair)}
  \item{discordance.by.subject}{a list of matrices (one for each subject) with the pair-wise discordance between the different genotyping instances of the subject}
  \item{correlation.by.subject}{a list of matrices (one for each
    subject) with the pair-wise correlation between the different
    genotyping instances of thesubject}

  If \code{corr.by.snp=TRUE}, \code{discordance.by.snp} will also have a
  column "correlation" with the correlation between duplicate subjects.
  For this calculation, the first two samples per subject are selected.
}
\author{Tushar Bhangale, Cathy Laurie, Stephanie Gogarten}

\seealso{\code{\link{GenotypeData}},
  \code{\link{duplicateDiscordanceAcrossDatasets}},
  \code{\link{duplicateDiscordanceProbability}},
  \code{\link{alleleFrequency}}
}
  
\examples{
library(GWASdata)
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)
data(illuminaScanADF)
genoData <-  GenotypeData(gds, scanAnnot=illuminaScanADF)

disc <- duplicateDiscordance(genoData, subjName.col="subjectID")

# minor allele discordance
afreq <- alleleFrequency(genoData)
minor.disc <- duplicateDiscordance(genoData, subjName.col="subjectID",
  minor.allele.only=TRUE, allele.freq=afreq[,"all"])

close(genoData)
}

\keyword{manip}
