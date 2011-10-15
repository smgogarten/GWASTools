\name{assocTestCPH}
\alias{assocTestCPH}
\title{Cox proportional hazards}
\description{Fits Cox proportional hazards model}

\usage{
assocTestCPH(genoData, event, time.to.event,
             covars, factor.covars,
             scan.chromosome.filter = NULL,
             scan.exclude = NULL,
             maf.filter = FALSE, 
             GxE = NULL, stratum = NULL,
             chromosome.set = NULL, block.size = 5000,
             verbose = TRUE,
             outfile = NULL)
}

\arguments{
  \item{genoData}{\code{\link{GenotypeData}} object, should contain sex and
    phenotypes in scan annotation}
  \item{event}{ name of scan variable in genoData for event to analyze}
  \item{time.to.event}{ name of scan variable in genoData for time to event}
  \item{covars}{vector of covariate terms for model
    (can include interactions as 'a:b', main effects correspond to
    scan variable names in genoData)}
  \item{factor.covars}{vector of names of covariates to be converted to
    factor}
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
  
  \item{maf.filter}{whether to filter results returned using \code{MAF*(1-MAF)
  > 75/(2*n)} where \code{MAF} = minor allele frequency and \code{n} = number of events}
  \item{GxE}{name of the covariate to use for E if genotype-by-environment (i.e. SNP:E) model is to be analyzed, in addition to the main effects (E can be a covariate interaction)}
  \item{stratum}{ name of variable to stratify on for a
    stratified analysis (use \code{NULL} if no stratified analysis
    needed)  }
  \item{chromosome.set}{integer vector with chromosome(s) to be analyzed.  Use 23, 24, 25, 26, 27 for X, XY, Y, M, Unmapped respectively.}
  \item{block.size}{ number of SNPs from a given chromosome to read in one block from genoData}
  \item{verbose}{Logical value specifying whether to show progress information.}
\item{outfile}{a character string to append in front of ".chr.i_k.RData" for naming the output data-frames; 
where i is the first chromosome, and k is the last chromosome used in
that call to the function.  "chr.i_k." will be
    omitted if \code{chromosome.set=NULL}.}
}
      
\details{
  This function performs Cox proportional hazards regression of a
  survival object (using the \code{\link{Surv}} function) on SNP
  genotype and other covariates.  It uses the \code{\link{coxph}}
  function from the R \code{\link{survival}} library.
  
  Individual samples can be included or excluded from the
  analysis using the \code{scan.exclude} parameter. Individual
  chromosomes can be included or excluded by specifying the indices of
  the chromosomes to be included in the \code{chromosome.set}
  parameter. Specific chromosomes for 
  specific samples can be included or excluded using the \code{scan.chromosome.filter}
    parameter.
  
  Both \code{scan.chromosome.filter} and \code{scan.exclude} may be used
  together. If a scan is excluded in EITHER, then it will be excluded
  from the analysis, but it does NOT need to be excluded in both. This
  design allows for easy filtering of anomalous scan-chromosome pairs
  using the \code{scan.chromosome.filter} matrix, but still allows easy
  exclusion of a specific group of scans (e.g. males or Caucasians)
  using \code{scan.exclude}.

  The argument \code{maf.filter} indicates whether to filter results returned
  using \code{2 * MAF * (1-MAF) * n > 75} where \code{MAF} = minor allele frequency and \code{n} = number of events. This filter was
  suggested by Ken Rice and Thomas Lumley, who found that without this requirement, at
  threshold levels of significance for genome-wide studies, Cox
  regression p-values based on standard asymptotic approximations can be
  notably anti-conservative.
}
      
\value{
  If \code{outfile=NULL} (default), all results are returned as a
  data.frame.  If \code{outfile} is specified, no data is returned
  but the function saves a data.frame with the naming convention as described by the argument \code{outfile}.
  Columns for the main effects model are:
    \item{index}{snp index}
    \item{snpID}{unique integer ID for SNP}
    \item{chr}{chromosome}
    \item{maf}{minor allele frequency calculated as appropriate for autosomal loci}
    \item{mafx}{minor allele frequency calculated as appropriate for X-linked loci}
    \item{beta}{ regression coefficient returned by the \code{\link{coxph}} function}
    \item{se}{ standard error of the regression coefficient returned by the \code{\link{coxph}} function}
    \item{z}{ z statistic returned by the \code{\link{coxph}} function}
    \item{pval}{ p-value for the z-statistic returned by the \code{\link{coxph}}
    function}
    \item{warned}{\code{TRUE} if a warning was issued}
    \item{n.events}{ number of events in complete cases for the given
      SNP}
    
  If \code{GxE} is not \code{NULL}, another data.frame is returned with the
  results of the genotype-by-environment model.  If \code{outfile=NULL},
  the function returns a list with names (\code{main},
  \code{GxE}); otherwise the GxE data.frame is saved as a separate output
  file.  Columns are:
    \item{index}{snp index}
    \item{snpID}{unique integer ID for SNP}
    \item{chr}{chromosome}
    \item{maf}{minor allele frequency calculated as appropriate for autosomal loci}
    \item{mafx}{minor allele frequency calculated as appropriate for X-linked loci}
    \item{warned}{\code{TRUE} if a warning was issued}
    \item{n.events}{ number of events in complete cases for the given SNP}
    \item{ge.lrtest}{Likelihood ratio test statistic for the GxE interaction}
    \item{ge.pval}{p-value for the likelihood ratio test statistic}
  
Warnings:

Another file will be saved with the name
"outfile.chr.i_k.warnings.RData" that contains any warnings generated by
the function.
}

\author{Cathy Laurie}

\seealso{\code{\link{GenotypeData}}, \code{\link{coxph}}}

\examples{
# an example of a scan chromosome matrix
# desiged to eliminate duplicated individuals
# and scans with missing values of sex
library(GWASdata)
data(affy_scan_annot)
scanAnnot <- ScanAnnotationDataFrame(affy_scan_annot)
samp.chr.matrix <- matrix(TRUE,nrow(scanAnnot),26)
dup <- duplicated(scanAnnot$subjectID)
samp.chr.matrix[dup | is.na(scanAnnot$sex),] <- FALSE
samp.chr.matrix[scanAnnot$sex=="F", 25] <- FALSE

# additionally, exclude YRI subjects
scan.exclude <- scanAnnot$scanID[scanAnnot$race == "YRI"]

# create some variables for the scans
scanAnnot$age <- rnorm(nrow(scanAnnot),mean=40, sd=10)
scanAnnot$event <- rbinom(nrow(scanAnnot),1,0.4)
scanAnnot$ttoe <- rnorm(nrow(scanAnnot),mean=100,sd=10)

# create data object
ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(ncfile)
genoData <-  GenotypeData(nc, scanAnnot=scanAnnot)

# variables
event <- "event"
time.to.event <- "ttoe"
covars <- c("sex", "age")
factor.covars <- "sex"

chr.set <- 21

res <- assocTestCPH(genoData,
  event="event", time.to.event="ttoe",
  covars=c("sex", "age"), factor.covars="sex",
  scan.chromosome.filter=samp.chr.matrix,
  scan.exclude=scan.exclude,
  chromosome.set=chr.set)

close(genoData)
}

\keyword{survival}
