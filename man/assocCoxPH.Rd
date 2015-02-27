\name{assocCoxPH}
\alias{assocCoxPH}
\title{Cox proportional hazards}
\description{Fits Cox proportional hazards model}

\usage{
assocCoxPH(genoData,
           event,
           time.to.event,
           gene.action = c("additive", "dominant", "recessive"),
           covar = NULL,
           ivar = NULL,
           strata = NULL,
           scan.exclude = NULL,
           effectAllele = c("minor", "alleleA"),
           snpStart = NULL,
           snpEnd = NULL,
           block.size = 5000,
           verbose = TRUE)
}

\arguments{
  \item{genoData}{a \code{\link{GenotypeData}} object}
  \item{event}{ name of scan annotation variable in \code{genoData} for event to analyze (chould be coded 0/1 or FALSE/TRUE)}
  \item{time.to.event}{ name of scan annotation variable in \code{genoData} for time to event}
  \item{gene.action}{
  "additive" coding sets the marker variable for homozygous minor allele samples = 2, heterozygous samples = 1, 
  and homozygous major allele samples = 0.
  "dominant" coding sets the marker variable for homozygous minor allele samples = 2, heterozygous samples = 2, 
  and homozygous major allele samples = 0.
  "recessive" coding sets the marker variable for homozygous minor allele samples = 2, heterozygous samples = 0, 
  and homozygous major allele samples = 0.
  (If \code{effectAllele="alleleA"}, the coding reflects alleleA instead of the minor allele.)
  }
  \item{covar}{
    a vector of the names of the covariates to adjust for (columns in the scan annotation of \code{genoData})
  }
  \item{ivar}{the name of the variable in \code{covar} to include as an interaction with genotype}
  \item{strata}{a vector of names of variables to stratify on for a
    stratified analysis}
  \item{scan.exclude}{a vector of scanIDs for scans to exclude}

  \item{effectAllele}{whether the effects should be returned in terms of the minor allele for the tested sample or the allele returned by \code{getAlleleA(genoData)}. If the minor allele is alleleB for a given SNP, this amounts to a sign change for the beta estimate.}

  \item{snpStart}{index of the first SNP to analyze, defaults to first SNP}
  \item{snpEnd}{index of the last SNP to analyze, defaults to last SNP}
  \item{block.size}{
    number of SNPs to read in at once.
  }
  \item{verbose}{logical for whether to print status updates}
}

\details{ 
  This function performs Cox proportional hazards regression of a
  survival object (using the \code{\link{Surv}} function) on SNP
  genotype and other covariates.  It uses the \code{\link{coxph}}
  function from the R \code{\link{survival}} library.
  
  It is recommended to filter results returned
  using \code{2 * MAF * (1-MAF) * n.events > 75} where \code{MAF} = minor allele frequency and \code{n.events} = number of events. This filter was
  suggested by Ken Rice and Thomas Lumley, who found that without this requirement, at
  threshold levels of significance for genome-wide studies, Cox
  regression p-values based on standard asymptotic approximations can be
  notably anti-conservative.

  Note: Y chromosome SNPs must be analyzed separately because they only use males.
}

\value{
a data.frame with the following columns
  \item{snpID}{the snpIDs}
  \item{chr}{chromosome SNPs are on}
  \item{n.events}{number of events in complete cases for each SNP}
  \item{MAF}{minor allele frequency}
  \item{minor.allele}{which allele is the minor allele }
  \item{Est}{beta estimate for genotype}
  \item{SE}{standard error of beta estimate for the genotype}
  \item{z.Stat}{z statistic for association}
  \item{z.pval}{p-value for association}

  if \code{ivar} is not \code{NULL}:
  \item{GxE.Stat}{Likelihood ratio test statistic for the genotype*ivar interaction parameter}
  \item{GxE.pval}{p-value for the likelihood ratio test statistic}
}

\author{Cathy Laurie, Matthew Conomos, Stephanie Gogarten}

\seealso{\code{\link{GenotypeData}}, \code{\link{coxph}}}

\examples{
library(GWASdata)
data(illuminaScanADF)
scanAnnot <- illuminaScanADF

# exclude duplicated subjects
scan.exclude <- scanAnnot$scanID[scanAnnot$duplicated]

# create some variables for the scans
scanAnnot$sex <- as.factor(scanAnnot$sex)
scanAnnot$age <- rnorm(nrow(scanAnnot), mean=40, sd=10)
scanAnnot$event <- rbinom(nrow(scanAnnot), 1, 0.4)
scanAnnot$ttoe <- rnorm(nrow(scanAnnot), mean=100, sd=10)

# create data object
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(gdsfile)
genoData <-  GenotypeData(gds, scanAnnot=scanAnnot)

res <- assocCoxPH(genoData,
                  event="event", time.to.event="ttoe",
  		  covar=c("sex", "age"),
  		  scan.exclude=scan.exclude,
  		  snpStart=1, snpEnd=100)

close(genoData)
}

\keyword{survival}
