\name{exactHWE}
\alias{exactHWE}
\title{Hardy-Weinberg Equilibrium testing}

\description{
This function performs exact Hardy-Weinberg Equilibrium testing (using Fisher's Test) over a selection of SNPs.  
It also counts genotype, calculates allele frequencies, and calculates inbreeding coefficients.
}

\usage{
exactHWE(genoData, 		
         scan.exclude = NULL,
         geno.counts = TRUE,
         snpStart = NULL,
         snpEnd = NULL,
         block.size = 5000, 
         verbose = TRUE)
}

\arguments{
  \item{genoData}{a \code{\link{GenotypeData}} object}
  
  \item{scan.exclude}{a vector of scanIDs for scans to exclude}
    
\item{geno.counts}{if \code{TRUE} (default), genotype counts are
  returned in the output data.frame.}

  \item{snpStart}{index of the first SNP to analyze, defaults to first SNP}
  \item{snpEnd}{index of the last SNP to analyze, defaults to last SNP}
  \item{block.size}{number of SNPs to read in at once}
  \item{verbose}{logical for whether to print status updates}
}

\details{
HWE calculations are performed with the \code{\link{HWExact}}
function in the \pkg{\link{GWASExactHW}} package.
  
For the X chromosome, only female samples will be used in all calculations (since males are excluded from HWE testing on this chromosome). 
The X chromosome may not be included in a block with SNPs from other chromosomes.
If the SNP selection includes the X chromosome, the scan annotation of genoData should include a "sex" column.  

Y and M and chromsome SNPs are not permitted in the SNP selection, since the HWE test is not valid for these SNPs.
}

\value{

a data.frame with the following columns
  \item{snpID}{the snpIDs}
  \item{chr}{chromosome SNPs are on}

If \code{geno.counts=TRUE}:
\item{nAA}{number of AA genotypes in samples}
\item{nAB}{number of AB genotypes in samples}
\item{nBB}{number of BB genotypes in samples} 

 \item{MAF}{minor allele frequency}
 \item{minor.allele}{which allele ("A" or "B") is the minor allele }
 \item{f}{the inbreeding coefficient}

 \item{pval}{exact Hardy-Weinberg Equilibrium (using Fisher's Test)
  p-value.  \code{pval} will be \code{NA} for monomorphic SNPs
  (\code{MAF=0}).}
}

\author{Ian Painter, Matthew P. Conomos, Stephanie Gogarten}

\seealso{\code{\link{HWExact}}}

\examples{
library(GWASdata)
data(illuminaScanADF)

# run only on YRI subjects
scan.exclude <- illuminaScanADF$scanID[illuminaScanADF$race != "YRI"]

# create data object
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(gdsfile)
genoData <-  GenotypeData(gds, scanAnnot=illuminaScanADF)
chr <- getChromosome(genoData)

# autosomal SNPs
auto <- range(which(is.element(chr, 1:22)))
hwe <- exactHWE(genoData, scan.exclude=scan.exclude, 
                snpStart=auto[1], snpEnd=auto[2])

# X chromosome SNPs must be run separately since they only use females
Xchr <- range(which(chr == 23))
hweX <- exactHWE(genoData, scan.exclude=scan.exclude, 
                 snpStart=Xchr[1], snpEnd=Xchr[2])

close(genoData)
}

\keyword{manip}
