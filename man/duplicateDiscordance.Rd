\name{duplicateDiscordance}
\alias{duplicateDiscordance}
\title{Duplicate discordance}
\description{
A function to compute all pair-wise genotype discordances between multiple genotyping instances of the same subject.
}
\usage{
duplicateDiscordance(genoData, subjName.col, 
                     scan.exclude = NULL, snp.exclude = NULL,
                     verbose = TRUE)
}
\arguments{
\item{genoData}{\code{\link{GenotypeData}} object} 
  \item{subjName.col}{A character string indicating the name of the annotation variable
  that will be identical for duplicate scans.}
\item{scan.exclude}{An integer vector containing the ids of scans to be
  excluded. }
\item{snp.exclude}{An integer vector containing the ids of SNPs to be excluded.}
\item{verbose}{Logical value specifying whether to show progress information.} 
}
\value{
  A list with the following components:
  \item{discordance.by.snp}{data frame with 5  columns: 1. snpID, 2. discordant (number of discordant pairs), 3. npair (number of pairs examined), 4. n.disc.subj (number of subjects with at least one discordance), 5. discord.rate (discordance rate i.e. discordant/npair)}
  \item{discordance.by.subject}{a list of matrices (one for each subject) with the pair-wise discordance between the different genotyping instances of the subject}
  \item{correlation.by.subject}{a list of matrices (one for each subject) with the pair-wise correlation between the different genotyping instances of the subject}
}
\author{Tushar Bhangale, Cathy Laurie}

\seealso{\code{\link{GenotypeData}},
  \code{\link{duplicateDiscordanceAcrossDatasets}},
  \code{\link{duplicateDiscordanceProbability}}
}
  
\examples{
library(GWASdata)
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(file)
data(affyScanADF)
genoData <-  GenotypeData(nc, scanAnnot=affyScanADF)

disc <- duplicateDiscordance(genoData, subjName.col="subjectID")
close(genoData)
}

\keyword{manip}
