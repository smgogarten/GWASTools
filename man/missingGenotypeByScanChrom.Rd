\name{missingGenotypeByScanChrom}
\alias{missingGenotypeByScanChrom}
\title{Missing Counts per Scan per Chromosome }
\description{
This function tabulates missing genotype calls for each scan for each chromosome.  
}
\usage{
missingGenotypeByScanChrom(genoData, snp.exclude = NULL, 
                           verbose = TRUE)
}


\arguments{
  \item{genoData}{\code{\link{GenotypeData}} object}
  \item{snp.exclude}{ A vector of IDs corresponding to the SNPs that should be excluded from the overall missing count.}
  \item{verbose}{Logical value specifying whether to show progress information.}
}



\details{
  This function calculates the percent of missing genotypes in each chromosome of each scan given in \code{genoData}.
  A "sex" variable must be present in the scan
  annotation slot of \code{genoData}.
}
 
 
\value{
This function returns a list with three components: "missing.counts,"
"snps.per.chr", and "missing.fraction."
  \item{missing.counts}{A matrix with rows corresponding to the scans and columns indicating unique chromosomes containing the number of missing SNP's for each scan and chromosome.}
  \item{snps.per.chr}{A vector containing the number of non-excluded
    SNPs for each chromosome.}
  \item{missing.fraction}{A vector containing the fraction of missing
    counts for each scan over all chromosomes, excluding the Y chromosome for females.}
}
\author{ Cathy Laurie }

\seealso{\code{\link{GenotypeData}}, \code{\link{missingGenotypeBySnpSex}}}


\examples{
library(GWASdata)
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(file)

# need scan annotation with sex
data(affy_scan_annot)
scanAnnot <- ScanAnnotationDataFrame(affy_scan_annot)
genoData <-  GenotypeData(nc, scanAnnot=scanAnnot)

missingRate <- missingGenotypeByScanChrom(genoData)
close(genoData)
}


\keyword{univar}