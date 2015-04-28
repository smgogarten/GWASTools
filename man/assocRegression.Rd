\name{assocRegression}
\alias{assocRegression}
\title{Association testing with regression}
\description{
  Run association testing with regression
}	

\usage{
assocRegression(genoData,
                outcome,
                model.type = c("linear", "logistic", "poisson", "firth"),
                gene.action = c("additive", "dominant", "recessive"),
                covar = NULL,
                ivar = NULL,
                scan.exclude = NULL,
	       	CI = 0.95,
                robust = FALSE,
                LRtest = FALSE,
                PPLtest = TRUE,
                effectAllele = c("minor", "alleleA"),
                snpStart = NULL,
                snpEnd = NULL,
                block.size = 5000,
                verbose = TRUE)
}

\arguments{
  \item{genoData}{a \code{\link{GenotypeData}} object}
  \item{outcome}{the name of the phenotype of interest (a column in the scan annotation of \code{genoData})}
  \item{model.type}{the type of model to be run.  "linear" uses \code{\link{lm}}, "logistic" uses \code{\link{glm}} with \code{family=binomial()}, "poisson" uses \code{\link{glm}} with \code{family=poisson()}, and "firth" uses \code{\link{logistf}}.}
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
  \item{scan.exclude}{a vector of scanIDs for scans to exclude}

  \item{CI}{a value between 0 and 1 defining the confidence level for the confidence interval calculations}

  \item{robust}{logical for whether to use sandwich-based robust standard errors for the "linear" or "logistic" method.  The default value is \code{FALSE}, and uses
  model based standard errors.  The standard error estimates are returned and also used for Wald Tests of significance.}

  \item{LRtest}{logical for whether to perform Likelihood Ratio Tests in addition to
  Wald tests (which are always performed).  
  Applies to linear, logistic, or poisson main effects only.
  NOTE:  Performing the LR tests adds a noticeable amount of computation time.}

  \item{PPLtest}{logical for whether to use the profile penalized likelihood to compute p values for the "firth" method (in addition to Wald tests, which are always performed).}

  \item{effectAllele}{whether the effects should be returned in terms of the minor allele for the tested sample (\code{effectAllele="minor"}) or the allele returned by \code{getAlleleA(genoData)}  (\code{effectAllele="alleleA"}). If the minor allele is alleleB for a given SNP, the difference between these two options will be a sign change for the beta estimate.}

  \item{snpStart}{index of the first SNP to analyze, defaults to first SNP}
  \item{snpEnd}{index of the last SNP to analyze, defaults to last SNP}
  \item{block.size}{number of SNPs to read in at once}
  \item{verbose}{logical for whether to print status updates}
}

\details{ 
  When using models without interaction terms, the association tests compare the model including the covariates and genotype value to the model 
  including only the covariates (a test of genotype effect).  When using a model with an interaction term, tests are performed for the 
  interaction term separately as well as a joint test of all the genotype terms (main effects and interactions) to detect any genotype effect.  
  All tests and p-values are always computed using Wald tests with p-values computed from Chi-Squared distribtuions.  
  The option of using either sandwich based robust standard errors (which make no model assumptions) or using model based standard errors for the 
  confidence intervals and Wald tests is specified by the \code{robust} parameter.  
  The option of also performing equivalent Likelihood Ratio tests is available and is specified by the \code{LRtest} parameter.
  
  For logistic regression models, if the SNP is monomorphic in either cases or controls, then the slope parameter is not well-defined, and the result will be \code{NA}.

  Note: Y chromosome SNPs must be analyzed separately because they only use males.
}

\value{
a data.frame with some or all of the following columns:
  \item{snpID}{the snpIDs}
  \item{chr}{chromosome SNPs are on}
  \item{effect.allele}{which allele ("A" or "B") is the effect allele}
  \item{EAF}{effect allele frequency}
  \item{MAF}{minor allele frequency}
  \item{n}{number of samples used to analyze each SNP}
  \item{n0}{number of controls (outcome=0) used to analyze each SNP}
  \item{n1}{number of cases (outcome=1) used to analyze each SNP}
  \item{Est}{beta estimate for genotype}
  \item{SE}{standard error of beta estimate for the genotype}
  \item{Wald.Stat}{chi-squared test statistic for association}
  \item{Wald.pval}{p-value for association}
  \item{LR.Stat}{likelihood ratio test statistic for association}
  \item{LR.pval}{p-value for association}
  \item{PPL.Stat}{profile penalized likelihood test statistic for association}
  \item{PPL.pval}{p-value for association}
  \item{GxE.Stat}{Wald test statistic for the genotype*ivar interaction parameter}
  \item{GxE.pval}{Wald test p-value for the genotype*ivar interaction parameter}
  \item{Joint.Stat}{Wald test statistic for jointly testing all genotype parameters}
  \item{Joint.pval}{Wald test p-value for jointly testing all genotype parameters}
}

\author{Tushar Bhangale, Matthew Conomos, Stephanie Gogarten}

\seealso{\code{\link{GenotypeData}}, \code{\link{lm}},
  \code{\link{glm}}, \code{\link{logistf}}, \code{\link{vcovHC}}, \code{\link{lrtest}}}

\examples{
library(GWASdata)
data(illuminaScanADF)
scanAnnot <- illuminaScanADF

# exclude duplicated subjects
scan.exclude <- scanAnnot$scanID[scanAnnot$duplicated]

# create some variables for the scans
scanAnnot$sex <- as.factor(scanAnnot$sex)
scanAnnot$age <- rnorm(nrow(scanAnnot), mean=40, sd=10)
scanAnnot$case.cntl.status <- rbinom(nrow(scanAnnot), 1, 0.4)
scanAnnot$blood.pressure[scanAnnot$case.cntl.status==1] <- rnorm(sum(scanAnnot$case.cntl.status==1), mean=100, sd=10)
scanAnnot$blood.pressure[scanAnnot$case.cntl.status==0] <- rnorm(sum(scanAnnot$case.cntl.status==0), mean=90, sd=5)

# create data object
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(gdsfile)
genoData <-  GenotypeData(gds, scanAnnot=scanAnnot)

## linear regression
res <- assocRegression(genoData,
		       outcome="blood.pressure",
                       model.type="linear",
                       covar=c("sex", "age"),
                       scan.exclude=scan.exclude,
 		       snpStart=1, snpEnd=100)

## logistic regression
res <- assocRegression(genoData,
		       outcome="case.cntl.status",
                       model.type="logistic",
                       covar=c("sex", "age"),
                       scan.exclude=scan.exclude,
 		       snpStart=1, snpEnd=100)

close(genoData)
}

\keyword{models}
\keyword{regression}


