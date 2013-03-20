\name{assocTestRegression}
\alias{assocTestRegression}
\title{Association tests}
\description{
This function performs regression based and likelihood ratio based association tests for both genotype main effects as well as interaction effects.  
It also computes genotype counts for association tests.
}
\usage{
assocTestRegression(genoData, outcome, model.type, 
                    covar.list = NULL, ivar.list = NULL,
                    gene.action.list = NULL, dosage = FALSE,
                    scan.chromosome.filter = NULL, 
                    scan.exclude = NULL, CI = 0.95,
                    robust = FALSE, LRtest = TRUE,
                    chromosome.set = NULL, block.set = NULL, 
                    block.size = 5000, verbose = TRUE,
                    outfile = NULL)
}

\arguments{
  \item{genoData}{\code{\link{GenotypeData}} object, should contain phenotypes and covariates in scan annotation}
  
  \item{outcome}{Vector (of length equal to the number of models) of
  names of the outcome variables for each model. These names must be in
  the scan annotation of \code{genoData}. e.g. c("case.cntl.status",
  "blood.pressure") will use "case.cntl.status" as the outcome for the
  first model and "blood pressure" for the second. Outcome variables
  must be coded as 0/1 for logistic regression.} 

  \item{model.type}{vector (of length equal to the number of models) with the types of models to be fitted.  
  The elements should be one of: "logistic", "linear", or "poisson".  e.g. c("logistic", "linear") will perform two tests: 
  the first using logistic regression, and the second using linear regression.}
  
  \item{covar.list}{list (of length equal to the number of models) of vectors containing the names of covariates to be used 
  in the regression model (blank, i.e. "" if none).  The default value is \code{NULL} and will include no covariates in any of the models.
  The covariate names must be in the scan annotation of \code{genoData}.
  e.g.
  \code{
  covar.list() <- list();
  covar.list[[1]] <- c("age","sex");
  covar.list[[2]] <- c("");
  }
  will use both "age" and "sex" as covariates for the first model and no covariates for the second model (this regresses on only the genotype).}

  \item{ivar.list}{list (of length equal to the number of models) of vectors containing the names of covariates for which to include an 
  interaction with genotype (blank, i.e. "" if none).  The default value is \code{NULL} and will include no interactions in any of the models.
  The covariate names must be in the scan annotation of \code{genoData}.
  e.g.
  \code{
  ivar.list() <- list();
  ivar.list[[1]] <- c("sex");
  ivar.list[[2]] <- c("");
  }
  will include a genotype*"sex" interaction term for the first model and no interactions for the second model.}

  \item{gene.action.list}{a list (of length equal to the number of models) of vectors containing the types of gene action models 
  to be used in the corresponding regression model.  Valid options are "additive", "dominant", and "recessive", referring to how 
  the minor allele is treated, as well as "dominance".
  "additive" coding sets the marker variable for homozygous minor allele samples = 2, heterozygous samples = 1, 
  and homozygous major allele samples = 0.
  "dominant" coding sets the marker variable for homozygous minor allele samples = 2, heterozygous samples = 2, 
  and homozygous major allele samples = 0.
  "recessive" coding sets the marker variable for homozygous minor allele samples = 2, heterozygous samples = 0, 
  and homozygous major allele samples = 0.
  "dominance" coding sets the marker variable for homozygous minor allele samples = major allele frequency, heterozygous samples = 0, 
  and homozygous major allele samples = minor allele frequency. This coding eliminates the additive component of variance for the 
  marker variable, leaving only the dominance component of variance.
  The default value is \code{NULL}, which assumes only an "additive" gene action model for every test.
  e.g.
  \code{
  gene.action.list() <- list();
  gene.action.list[[1]] <- c("additive");
  gene.action.list[[2]] <- c("dominant", "recessive");
  }
  will run the first model using "additive" gene action, and will run the second model using both "dominant" and "recessive" gene actions.}

  \item{dosage}{logical for whether or not the genotype values are imputed dosages.  The default value is \code{FALSE} for true genotype calls.
  When using imputed dosages, the \code{gene.action} must be additive, and genotype counts will not be calculated.}
  
  \item{scan.chromosome.filter}{a logical matrix that can be used to
  exclude some chromosomes, some scans, or some specific scan-chromosome
  pairs. Entries should be \code{TRUE} if that scan-chromosome pair should be
  included in the analysis, \code{FALSE} if not. The number of rows must be
  equal to the number of scans in \code{genoData}, and the number of
  columns must be equal to the largest integer chromosome value in \code{genoData}. The column number must match the chromosome number.
e.g. A scan.chromosome.filter matrix used for an analyis when
  \code{genoData} has SNPs with chromosome=(1-24, 26, 27) (i.e. no Y
  (25) chromosome SNPs) must have 27 columns (all \code{FALSE} in the
  25th column). But a scan.chromosome.filter matrix used for an analysis
  \code{genoData} has SNPs chromosome=(1-26) (i.e no Unmapped (27)
  chromosome SNPs) must have only 26 columns.}

  \item{scan.exclude}{an integer vector containing the IDs of entire scans to
    be excluded.}
  
  \item{CI}{sets the confidence level for the confidence interval calculations.  Confidence intervals are computed at every SNP;
  for the odds ratio when using logistic regression, for the linear trend parameter when using linear regression, and for the rate ratio
  when using Poisson regression.  The default value is 0.95 (i.e. a 95\% confidence interval). The confidence level must be between 0 and 1.}

  \item{robust}{logical for whether to use sandwich-based robust standard errors.  The default value is \code{FALSE}, and uses
  model based standard errors.  The standard error estimates are returned and also used for Wald Tests of significance.}
  
  \item{LRtest}{logical for whether to perform Likelihood Ratio Tests.  The default value is \code{TRUE}, and performs LR tests in addition to
  Wald tests (which are always performed).  NOTE:  Performing the LR tests adds a noticeable amount of computation time.}

  \item{chromosome.set}{integer vector with chromosome(s) to be
  analyzed.  Use 23, 24, 25, 26, 27 for X, XY, Y, M, Unmapped
  respectively.}

  \item{block.set}{list (of length equal to \code{length(chromosome.set)}) of vectors where every vectors contains the indices of the 
  SNP blocks (on that chromosome) to be analyzed.
  e.g.
  \code{chromosome.set <- c(1,2); block.set <- list(); chr.1 <- c(1,2,3);
  chr.2 <- c(5,6,7,8); block.set$chr.1 <- chr.1; block.set$chr.2 <- chr.2;}
  will analyze first three block on chromosome 1 and 5th through 8th 
  blocks on chromosome 2. The actual number of SNPs analyzed will depend on \code{block.size}. 
  Default value is \code{NULL}. If \code{block.set == NULL}, all the SNPs on chromosomes in \code{chromosome.set} will be analyzed. }
  
  \item{block.size}{Number of SNPs to be read from \code{genoData} at once.}
  
  \item{verbose}{if \code{TRUE} (default), will print status updates while the function runs.
  e.g. it will print "chr 1 block 1 of 10" etc. in the R console after each block of SNPs is done being analyzed.}

  \item{outfile}{a character string to append in front of
  ".model.j.gene_action.chr.i_k.RData" for naming the output data-frames; where j is the model number, gene_action is the gene.action type, 
  i is the first chromosome, and k is the last chromosome used in that call to the function.  "chr.i_k." will be omitted if 
  \code{chromosome.set=NULL}.  If set to \code{NULL} (default), then the results are returned to the R console.}
}

\details{
  When using models without interaction terms, the association tests compare the model including the
  covariates and genotype value to the model including only the covariates (a test of genotype effect).  When using a model with interaction terms,
  the association tests compare the model including all the interactions, the covariates, and the genotype value to the model
  with only the covariates and genotype value (jointly testing for all the interaction effects).  All tests and p-values are always
  computed using Wald tests with p-values computed from Chi-Squared distribtuions.  The option of using either sandwich based robust standard errors (which make no model assumptions)
  or using model based standard errors for the confidence intervals and Wald tests is specified by the \code{robust} parameter.
  The option of also performing equivalent Likelihood Ratio tests is available and is specified by the \code{LRtest} parameter.
  
  Three types of regression models are available: "logistic", "linear", or "poisson".
  Multiple models can be run at the same time by putting multiple arguments in the \code{outcome}, \code{model.type},
  \code{covar.list}, \code{ivar.list}, and \code{gene.action.list} parameters.  For each model, available gene action models 
  are "additive", "dominant", "recessive", and "dominance."  See above for the correct usage of each of these.
  
  For logistic regression models, if the SNP is monomorphic in either cases or controls, then the slope parameter is not well-defined.
  In this situation, an error message will be returned (see \code{model.N.gene_action.warningOrError} in the Value section below for details), 
  and the regression of this SNP will not be performed.  If a test of significance is still desired for these SNPs, 
  we suggest performing either a Fisher's Exact Test using the \code{assocTestFisherExact} function provided in GWASTools 
  or performing a trend test (using \code{model.type = "linear"} in this function).

  Individual samples can be included or excluded from the
  analysis using the \code{scan.exclude} parameter. Individual
  chromosomes can be included or excluded by specifying the indices of
  the chromosomes to be included in the \code{chromosome.set}
  parameter.  Specific chromosomes for specific samples can be included or excluded using the \code{scan.chromosome.filter} parameter.  
  The inclusion or exclusion of specific blocks of SNP's on each chromosome can be specified using the \code{block.set}
  parameter. Note that the actual SNP's included or excluded will change according to the value of \code{block.size}.
  
  Both \code{scan.chromosome.filter} and \code{scan.exclude} may be used
  together. If a scan is excluded in EITHER, then it will be excluded
  from the analysis, but it does NOT need to be excluded in both. This
  design allows for easy filtering of anomalous scan-chromosome pairs
  using the \code{scan.chromosome.filter} matrix, but still allows easy
  exclusion of a specific group of scans (e.g. males or Caucasians)
  using \code{scan.exclude}.
  
  This function allows for the usage of imputed dosages in place of genotypes in the additive model.
}

\value{
  If \code{outfile=NULL} (default), all results are returned as a single
  data.frame.  If \code{outfile} is specified, no data is returned but
  the function saves a data-frame for each model gene-action pair, with
  the naming convention as described by the variable \code{outfile}.

  The first column of each data-frame is: 
  \item{snpID}{snpID (from \code{genoData}) of the SNP}
  
  After these first three columns, for every model gene-action pair there are the following columns:
  Here, "model.M" is the name assigned to the test where M = 1, 2,..., length(model.type), 
  and "gene_action" is the gene-action type of the test (one of "additive", "dominant", "recessive", or "dominance").
  \item{model.M.n}{sample size for the regression}
  
  For tests that use linear regression (will be \code{NA} if using imputed dosages for genotypes):
  \item{model.M.nAA}{ 	number of AA genotypes in samples}
  \item{model.M.nAB}{ 	number of AB genotypes in samples}
  \item{model.M.nBB}{ 	number of BB genotypes in samples}

  For tests that use logistic regression (will be \code{NA} if using imputed dosages for genotypes):
  \item{model.M.nAA.cc0}{number of AA genotypes in samples with outcome coded as 0}
  \item{model.M.nAB.cc0}{number of AB genotypes in samples with outcome coded as 0}
  \item{model.M.nBB.cc0}{number of BB genotypes in samples with outcome coded as 0}
  \item{model.M.nAA.cc1}{number of AA genotypes in samples with outcome coded as 1}
  \item{model.M.nAB.cc1}{number of AB genotypes in samples with outcome coded as 1}
  \item{model.M.nBB.cc1}{number of BB genotypes in samples with outcome coded as 1}
  
  \item{model.M.MAF}{minor allele frequency. Note that calculation of allele frequency for the X chromosome is different than that for the 
  autosomes and the XY (pseudo-autosomal) region. Hence if chromosome.set includes 23, \code{genoData} should provide the sex of the 
  scan ("M" or "F") i.e. there should be a column named "sex" with "F" for females and "M" for males.}
  \item{model.M.minor.allele}{the minor allele. Takes values "A" or "B".}
   
  \item{model.M.gene_action.warningOrError}{report of different possible warnings or errors:  
  0 if controls are monomorphic (logistic regression only), 1 if cases are monomorphic (logistic refression only), 2 if all samples are monomorphic, 
  9 if a warning or error occured during model fitting, \code{NA} if none}
  \item{model.M.gene_action.Est.G}{estimate of the regression coefficient for the genotype term. 
  See the description in \code{gene.action.list} above for interpretation.}
  \item{model.M.gene_action.SE.G}{standard error of the regression coefficient estimate for the genotype term.  
  Could be either sandwich based (robust) or model based; see description in \code{robust}.}
  
  For tests with interaction variables:
  Here, "ivar_name" refers to the name of the interaction variable; if there are multiple interaction variables, 
  there will be a column with each different "ivar_name".
  \item{model.M.gene_action.Est.G.ivar_name}{estimate of the regression coefficient for the interaction between genotype and ivar_name.}
  \item{model.M.gene_action.SE.G.ivar_name}{standard error of the regression coefficient estimate.
  Could be either sandwich based (robust) or model based; see description in \code{robust}.}
  
  For tests that use logistic regression with no interaction variables:
  \item{model.M.gene_action.OR.G}{odds ratio for the genotype term. This is exp(the regression coefficient). 
  See the description in "gene.action.list" above for interpretation.}
  \item{model.M.gene_action.OR_L95.G}{lower 95\% confidence limit for the odds ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  \item{model.M.gene_action.OR_U95.G}{upper 95\% confidence limit for the odds ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  
  For tests that use logistic regression and interaction variables:
  \item{model.M.gene_action.OR.G.ivar_name}{odds ratio for the genotype*ivar_name interaction term. 
  This is exp(the interaction regression coefficient).  A separate odds ratio is calculated for each interaction term.  
  See the description in "gene.action.list" above for interpretation.}
  \item{model.M.gene_action.OR_L95.G.ivar_name}{lower 95\% confidence limit for the odds ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  \item{model.M.gene_action.OR_U95.G.ivar_name}{upper 95\% confidence limit for the odds ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  
  For tests that use linear regression with no interaction variables:
  \item{model.M.gene_action.L95.G}{lower 95\% confidence limit for the genotype coefficient 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  \item{model.M.gene_action.U95.G}{upper 95\% confidence limit for the genotype coefficient
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  
  For tests that use linear regression and interaction variables:
  \item{model.M.gene_action.L95.G.ivar_name}{lower 95\% confidence limit for the genotype*ivar_name interaction coefficient 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  \item{model.M.gene_action.U95.G.ivar_name}{upper 95\% confidence limit for the genotype*ivar_name interaction coefficient
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  
  For tests that use Poisson regression with no interaction variables:
  \item{model.M.gene_action.RR.G}{rate ratio for the genotype term. This is exp(the regression coefficient). 
  See the description in "gene.action.list" above for interpretation.}
  \item{model.M.gene_action.RR_L95.G}{lower 95\% confidence limit for the rate ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  \item{model.M.gene_action.RR_U95.G}{upper 95\% confidence limit for the rate ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  
  For tests that use Poisson regression and interaction variables:
  \item{model.M.gene_action.RR.G.ivar_name}{rate ratio for the genotype*ivar_name interaction term. 
  This is exp(the interaction regression coefficient).  A separate odds ratio is calculated for each interaction term.  
  See the description in "gene.action.list" above for interpretation.}
  \item{model.M.gene_action.RR_L95.G.ivar_name}{lower 95\% confidence limit for the rate ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  \item{model.M.gene_action.RR_U95.G.ivar_name}{upper 95\% confidence limit for the rate ratio 
  (95 will be replaced with whatever confidence level is chosen in \code{CI}).}
  
  For tests with no interaction variables:
  \item{model.M.gene_action.Wald.Stat.G}{value of the Wald test statistic for testing the genotype parameter}
  \item{model.M.gene_action.Wald.pval.G}{Wald test p-value, calculated from a Chi-Squared distribution.  
  This can be calculated using either sandwich based robust standard errors or model based standard errors (see \code{robust}).}
  
  For tests with interaction variables:
  \item{model.M.gene_action.Wald.Stat.GxE}{value of the Wald test statistic for jointly testing all genotype interaction parameters}
  \item{model.M.gene_action.Wald.pval.GxE}{Wald test p-value for jointly testing all genotype interaction parameters, calculated from a Chi-Squared distribution.  
  This can be calculated using either sandwich based robust standard errors or model based standard errors (see \code{robust}).}
  
  If \code{LRtest = TRUE}, for tests with no interaction variables:
  \item{model.M.gene_action.LR.Stat.G}{value of the Likelihood Ratio test statistic for testing the genotype parameter}
  \item{model.M.gene_action.LR.pval.G}{Likelihood Ratio test p-value.}
  
  If \code{LRtest = TRUE}, for tests with interaction variables:
  \item{model.M.gene_action.LR.Stat.GxE}{value of the Likelihood Ratio test statistic for jointly testing all genotype interaction parameters}
  \item{model.M.gene_action.LR.pval.GxE}{Likelihood Ratio test p-value for jointly testing all genotype interaction parameters.}
  


Attributes:

There is also an attribute for each output data-frame called "model" that shows the model used for the test. 
This can be viewed with the following R command:
\code{attr(mod.res, "model")}
where mod.res is the output data-frame from the function.
The \code{attr()} command will return something like:
model.1.additive
"case.cntl.status ~ age + sex , logistic regression, additive gene action"

There is another attribute called "SE" that shows if Robust or Model Based standard errors were used for the test.
This can be viewed with the following R command:
\code{attr(mod.res, "SE")}
where mod.res is the output data-frame from the function.


Warnings:

Another file will be saved with the name "outfile.chr.i_k.warnings.RData" that contains any warnings generated by the function.
An example of what would be contained in this file:
Warning messages:
1: Model 1 , Y chromosome tests are confounded with sex and should be run separately without sex in the model
2: Model 2 , Y chromosome tests are confounded with sex and should be run separately without sex in the model 
}

\author{Matthew P. Conomos, Tushar Bhangale}

\seealso{\code{\link{GenotypeData}}, \code{\link{lm}},
  \code{\link{glm}}, \code{\link{vcov}}, \code{\link{vcovHC}}, \code{\link{lrtest}}}

\examples{
# The following example would perform 3 tests (from 2 models): 
# the first a logistic regression of case.cntl.status on genotype, age, and sex, including an interaction term between genotype and sex, using additive gene action; 
# the second a linear regression of blood pressure on genotype using dominant gene action, 
# and the third, a linear regression of blood pressure on genotype again, but this time using recessive gene action.
# This test would only use chromosome 21.  
# It would perform both robust Wald tests using sandwich based robust standard errors as well as Likelihood Ratio tests.

# an example of a scan chromosome matrix
# desiged to eliminate duplicated individuals
# and scans with missing values of sex
library(GWASdata)
data(affy_scan_annot)
scanAnnot <- ScanAnnotationDataFrame(affy_scan_annot)
samp.chr.matrix <- matrix(TRUE,nrow(scanAnnot),26)
dup <- duplicated(scanAnnot$subjectID)
samp.chr.matrix[dup | is.na(scanAnnot$sex),] <- FALSE

# additionally, exclude YRI subjects
scan.exclude <- scanAnnot$scanID[scanAnnot$race == "YRI"]

# create some variables for the scans
scanAnnot$sex <- as.factor(scanAnnot$sex)
scanAnnot$age <- rnorm(nrow(scanAnnot),mean=40, sd=10)
scanAnnot$case.cntl.status <- rbinom(nrow(scanAnnot),1,0.4)
scanAnnot$blood.pressure[scanAnnot$case.cntl.status==1] <- rnorm(sum(scanAnnot$case.cntl.status==1),mean=100,sd=10)
scanAnnot$blood.pressure[scanAnnot$case.cntl.status==0] <- rnorm(sum(scanAnnot$case.cntl.status==0),mean=90,sd=5)

# create data object
ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(ncfile)
genoData <-  GenotypeData(nc, scanAnnot=scanAnnot)

# set regression variables and models
outcome <- c("case.cntl.status","blood.pressure")

covar.list <- list()
covar.list[[1]] <- c("age","sex")
covar.list[[2]] <- c("")

ivar.list <- list();
ivar.list[[1]] <- c("sex");
ivar.list[[2]] <- c("");
  
model.type <- c("logistic","linear")

gene.action.list <- list()
gene.action.list[[1]] <- c("additive")
gene.action.list[[2]] <- c("dominant", "recessive")

chr.set <- 21

outfile <- tempfile()

assocTestRegression(genoData,
                    outcome = outcome,
					model.type = model.type,
                    covar.list = covar.list,
                    ivar.list = ivar.list,
                    gene.action.list = gene.action.list,
                    scan.chromosome.filter = samp.chr.matrix,
                    scan.exclude = scan.exclude,
                    CI = 0.95,
                    robust = TRUE,
					LRtest = TRUE,
                    chromosome.set = chr.set,
                    outfile = outfile)

model1 <- getobj(paste(outfile, ".model.1.additive.chr.21_21.RData", sep=""))
model2 <- getobj(paste(outfile, ".model.2.dominant.chr.21_21.RData", sep=""))
model3 <- getobj(paste(outfile, ".model.2.recessive.chr.21_21.RData", sep=""))

close(genoData)
unlink(paste(outfile, "*", sep=""))

# In order to run the test on all chromosomes, it is suggested to run the function in parallel.
# To run the function in parallel the following unix can be used:
# R --vanilla --args 21 22 < assoc.analysis.r >logfile.txt &
# where the file assoc.analysis.r will include commands similar to this example 
# where chromosome.set and/or block.set can be passed to R using --args
# Here, tests on chromosomes 21 and 22 are performed; these could be replaced by any set of chromosomes
# these values are retrieved in R by putting a
# chr.set <- as.numeric(commandArgs(trailingOnly=TRUE))
# command in assoc.analysis.r 
}

\keyword{models}
\keyword{regression}


