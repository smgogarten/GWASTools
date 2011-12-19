\name{hetBySnpSex}
\alias{hetBySnpSex}
\title{Heterozygosity by SNP and sex
}
\description{This function calculates the percent of
  heterozygous genotypes for males and females for each SNP.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
hetBySnpSex(genoData, scan.exclude = NULL, 
            verbose = TRUE)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{
\item{genoData}{\code{\link{GenotypeData}} object}

\item{scan.exclude}{ An integer vector containing the id's of scans to be excluded. }

\item{verbose}{Logical value specifying whether to show progress information.}
  
 
}





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
This function calculates the percent of heterozygous
genotypes for males and females for each SNP given in \code{genoData}. A
"sex"  variable must be present in the scan
  annotation slot of \code{genoData}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{
The result is a matrix containing the heterozygosity rates with snps as
rows and 2 columns ("M" for males and "F" for females).
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\author{Cathy Laurie
}
\seealso{\code{\link{GenotypeData}}, \code{\link{hetByScanChrom}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
library(GWASdata)
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(file)

# need scan annotation with sex
data(affyScanADF)
genoData <- GenotypeData(nc, scanAnnot=affyScanADF)

het <- hetBySnpSex(genoData)
close(genoData)
}
\keyword{ manip }


