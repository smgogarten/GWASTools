\name{duplicateDiscordanceAcrossDatasets}
\alias{duplicateDiscordanceAcrossDatasets}
\title{Duplicate discordance across datasets}

\description{
  Finds number of discordant genotypes by SNP in pairs of duplicate scans
  of the same subject across multiple datasets.  
}

\usage{
duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
  subjName.cols, snpName.cols,
  alleleA.cols=NULL, alleleB.cols=NULL,
  one.pair.per.subj=TRUE, minor.allele.only=FALSE, 
  scan.exclude1=NULL, scan.exclude2=NULL, snp.include=NULL,
  verbose=TRUE)
}

\arguments{
  \item{genoData1}{\code{\link{GenotypeData}} object containing the
    first dataset.}
  \item{genoData2}{\code{\link{GenotypeData}} object containing the
    second dataset.}
  \item{subjName.cols}{2-element character vector indicating the names of the
    annotation variables that will be identical for duplicate scans in
    the two datasets.}
  \item{snpName.cols}{2-element character vector indicating the names of the
    annotation variables that will be identical for the same SNPs in the
    two datasets.}
  \item{alleleA.cols}{2-element character vectors indicating the names
  of the annotation variables with allele A for the two datasets.}
  \item{alleleB.cols}{2-element character vectors indicating the names
  of the annotation variables with allele B for the two datasets.}
\item{one.pair.per.subj}{A logical indicating whether a single pair of
  scans should be randomly selected for each subject with more than 2 scans.}
\item{minor.allele.only}{A logical indicating whether discordance should be calculated
  only between pairs of scans in which at least one scan has a genotype
  with the minor allele (i.e., exclude major allele homozygotes).}
\item{scan.exclude1}{An integer vector containing the ids of scans to be
  excluded from the first dataset. }
\item{scan.exclude2}{An integer vector containing the ids of scans to be
  excluded from the second dataset. }
  \item{snp.include}{List of SNPs to include in the comparison.  Should
    match the contents of the columns referred to by \code{snpName.cols}.}
\item{verbose}{Logical value specifying whether to show progress information.} 
}

\details{
  \code{duplicateDiscordanceAcrossDatasets} calculates discordance metrics both by
  scan and by SNP.  If \code{one.pair.per.subj=TRUE} (the default), each
  subject with more than two duplicate genotyping instances will have
  one scan from each dataset randomly selected for computing discordance.  If
  \code{one.pair.per.subj=FALSE}, discordances will be calculated
  pair-wise for all possible cross-dataset pairs for each subject.

  If \code{minor.allele.only=TRUE}, the allele frequency will be
  calculated in \code{genoData1}, using only samples common to both datasets.
  
  If \code{snp.include=NULL} (the default), discordances will be found
  for all SNPs common to both datasets.
}

\value{
  A list with the following components:
  \item{discordance.by.snp}{data frame with 4 columns: 1. discordant
  (number of discordant pairs), 2. npair (number of pairs examined),
  3. n.disc.subj (number of subjects with at least one discordance),
  4. discord.rate (discordance rate i.e. discordant/npair).  Row names
  are the common snp ID.}
  \item{discordance.by.subject}{a list of matrices (one for each subject) with the pair-wise discordance between the different genotyping instances of the subject}
  
  If no duplicate scans or no common SNPs are found, issues a warning
  message and returns \code{NULL}.
}

\author{Stephanie Gogarten, Jess Shen}

\seealso{\code{\link{GenotypeData}}, \code{\link{duplicateDiscordance}},
  \code{\link{duplicateDiscordanceProbability}}}

\examples{
library(GWASdata)

# dataset 1
file1 <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc1 <- NcdfGenotypeReader(file1)
data(affySnpADF)
data(affyScanADF)
data1 <-  GenotypeData(nc1, snpAnnot=affySnpADF, scanAnnot=affyScanADF)

# dataset 2
file2 <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
nc2 <- NcdfGenotypeReader(file2)
data(illuminaSnpADF)
data(illuminaScanADF)
data2 <-  GenotypeData(nc2, snpAnnot=illuminaSnpADF, scanAnnot=illuminaScanADF)

discord <- duplicateDiscordanceAcrossDatasets(data1, data2,
             subjName.cols=c("CoriellID", "CoriellID"),
             snpName.cols=c("rsID", "rsID"))
close(data1)
close(data2)
}

\keyword{manip}
