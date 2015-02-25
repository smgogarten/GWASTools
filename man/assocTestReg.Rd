\name{assocTestReg}
\alias{assocTestReg}
\title{Association testing with regression}
\description{
  Run association testing with regression
}	

\usage{
assocTestReg(genoData,
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
  \item{outcome}{the name of the phenotype of interest}
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
    a vector of the names of the covariates to adjust for
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

  \item{effectAllele}{whether the effects should be returned in terms of the minor allele for the tested sample or the allele returned by \code{getAlleleA(genoData)}. If the minor allele is alleleB for a given SNP, this amounts to a sign change for the beta estimate.}

  \item{snpStart}{index of the first SNP to analyze, defaults to first SNP}
  \item{snpEnd}{index of the last SNP to analyze, defaults to last SNP}
  \item{block.size}{
    number of SNPs to read in at once.
  }
  \item{verbose}{logical for whether to print status updates}
}

\details{ 
  Note: Y chromosome SNPs must be analyzed separately b/c they only use males.
}

\value{
a data.frame with the following columns
  \item{snpID}{the snpIDs}
  \item{chr}{chromosome SNPs are on}
  \item{n}{number of samples used to analyze each SNP}
  \item{MAF}{minor allele frequency}
  \item{minor.allele}{which allele is the minor allele }
  \item{Est}{beta estimate for genotype}
  \item{SE}{standard error of beta estimate for the genotype}
  \item{Wald.Stat}{chi-squared test statistic for association}
  \item{Wald.pval}{p-value for association}

  if \code{model.type="linear" or "logistic"} and \code{LRtest=TRUE}:
  \item{LR.Stat}{likelihood ratio test statistic for association}
  \item{LR.pval}{p-value for association}

  if \code{model.type="firth"} and \code{PPLtest=TRUE}:
  \item{PPL.Stat}{profile penalized likelihood test statistic for association}
  \item{PPL.pval}{p-value for association}

  if \code{model.type="linear" or "logistic"} and \code{ivar} is not \code{NULL}:
  \item{GxE.Stat}{Wald test statistic for the genotype*ivar interaction parameter}
  \item{GxE.pval}{Wald test p-value for the genotype*ivar interaction parameter}
  \item{Joint.Stat}{Wald test statistic for jointly testing all genotype parameters}
  \item{Joint.pval}{Wald test p-value for jointly testing all genotype parameters}

}

\author{Matt Conomos, Stephanie Gogarten}
