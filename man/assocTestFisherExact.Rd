\name{assocTestFisherExact}
\alias{assocTestFisherExact}
\title{Association tests}
\description{
This function performs Fisher's Exact Test using allele counts for cases and controls.  
It takes the output from \code{\link{assocTestRegression}} as its input.
}
\usage{
assocTestFisherExact(dat, outfile = NULL)
}

\arguments{
  \item{dat}{a data.frame of output from \code{\link{assocTestRegression}} run with \code{model.type = "logistic"} (a case/control test).
  It should contain all columns of the output and only the rows (SNPs) that the user wishes to perform Fisher's Exact Test on.}
 
  \item{outfile}{a character string to append in front of "FisherExact.Rdata" for naming the output data-frame. 
  If set to \code{NULL} (default), then the results are returned to the R console.}
}

\details{
  This function performs a basic Fisher's Exact Test to test for differences in allele frequencies between cases and controls;
  it compares the "A" and "B" allele frequencies between cases and controls.
  
  This function uses the output from \code{\link{assocTestRegression}} run with \code{model.type = "logistic"} as its input.
  It uses the output genotype counts for cases and controls, converts them to allele counts and performs the Fisher's Exact Test
  to calculate an allelic odds ratio (the odds of being a case for the minor allele compared to the major allele), a 95\% confidence interval, 
  and a p-value.  
  
  One suggested use of this function is to perform significance tests on SNPs that are monomorphic in either cases or controls,
  as a standard logistic regression test is not well-defined in this case.  The \code{\link{assocTestRegression}} function will return an
  error for these SNPs; see its help page for more detail.
}

\value{
  If \code{outfile=NULL} (default), all results are returned as a data.frame.  
  If \code{outfile} is specified, no data is returned but the function saves 
  a data-frame with the naming convention as described by the variable \code{outfile}.

  The first five columns of the data-frame are taken from \code{dat}: 
  \item{snpID}{snpID of the SNP}
  \item{MAF}{minor allele frequency. Note that calculation of allele frequency for the X chromosome is different than that for the 
  autosomes and the XY (pseudo-autosomal) region.}
  \item{minor.allele}{the minor allele. Takes values "A" or "B".}
  \item{n}{sample size for the regression}
  \item{regression.warningOrError}{report of different possible warnings or errors from the regression test:  
  0 if controls are monomorphic (logistic regression only), 1 if cases are monomorphic (logistic refression only), 
  9 if a warning or error occured during model fitting, \code{NA} if none}
  
  \item{Fisher.OR}{odds ratio from the Fisher's Exact test of allele counts. It is the odds of being a case for the minor allele compared to the major allele.}
  \item{Fisher.OR_L95}{lower 95\% confidence limit for the odds ratio.}
  \item{Fisher.OR_U95}{upper 95\% confidence limit for the odds ratio.}
  \item{Fisher.pval}{Fisher's Exact test p-value.}
  
  \item{nA.cc0}{number of A alleles among samples with outcome coded as 0}
  \item{nB.cc0}{number of B alleles among samples with outcome coded as 0}
  \item{nA.cc1}{number of A alleles among samples with outcome coded as 1}
  \item{nB.cc1}{number of B alleles among samples with outcome coded as 1}
}
  
\author{Matthew P. Conomos}

\seealso{\code{\link{assocTestRegression}}}

\examples{
# The following example would take the output from association tests run on chromosome 22 using assocTestRegression
# and perform the Fisher's Exact Test on those that were monomorphic in either the cases or the controls.
# The output would be saved as "chr22test.FisherExact.RData"

# run assocTestRegression
library(GWASdata)
data(affyScanADF)
ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(ncfile)
genoData <-  GenotypeData(nc, scanAnnot=affyScanADF)
mydat <- assocTestRegression(genoData, outcome="status",
  model.type="logistic", chromosome.set=22)

# subset rows of those SNPs that are monomorphic in cases or controls; keep all columns
mono.dat <- mydat[which(mydat$model.1.additive.warningOrError == 0 |
                        mydat$model.1.additive.warningOrError ==1),]

# perform the Fisher's Exact Test
assocTestFisherExact(dat = mono.dat, outfile = "chr22test")

# load the output
outfile <- "chr22test.FisherExact.RData"
fisher.res <- getobj(outfile)
head(fisher.res)
unlink(outfile)
}

