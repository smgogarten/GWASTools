\name{hetByScanChrom}
\alias{hetByScanChrom}
\title{Heterozygosity rates by scan and chromosome
}
\description{This function calculates the fraction of heterozygous genotypes for each chromosome for a set of scans.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
hetByScanChrom(genoData, snp.exclude = NULL, 
               verbose = TRUE)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{
\item{genoData}{\code{\link{GenotypeData}} object.  Chromosomes are
  expected to be in contiguous blocks.}

\item{snp.exclude}{ An integer vector containing the id's of SNPs to be excluded. }

\item{verbose}{Logical value specifying whether to show progress information.}
  
 
}





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{

This function calculates the percent of heterozygous and missing genotypes in each chromosome of each scan given in \code{genoData}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{
The result is a matrix containing the heterozygosity rates with scans as
  rows and chromosomes as columns, including a column "A" for all autosomes.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\author{Cathy Laurie
}
\seealso{\code{\link{GenotypeData}}, \code{\link{hetBySnpSex}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)
genoData <- GenotypeData(gds)
het <- hetByScanChrom(genoData)
close(genoData)
}
\keyword{ manip }


