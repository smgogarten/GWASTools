\name{simulateGenotypeMatrix}
\alias{simulateGenotypeMatrix}
\title{ Simulate Genotype Matrix & Load into NetCDF File}
\description{
This function creates a netCDF file with dimensions 'snp' and 'sample' and variables 'sampleID', 'genotype', 'position' and 'chromosome'.  These variables hold simulated data as described below.  
Mainly, this function is intended to be used in examples involving genotype matrices.
}
\usage{
simulateGenotypeMatrix(n.snps=10, n.chromosomes=10,
                       n.samples=1000, ncdf.filename,
                       silent=TRUE)
}
\arguments{
  \item{n.snps}{ An integer corresponding to the number of SNPs per chromosome, the default value is 10.  For this function, the number of SNPs is assumed to be the same for every chromosome. }
  \item{n.chromosomes}{ An integer value describing the total number of chromosomes with default value 10. }
  \item{n.samples}{ An integer representing the number of samples for our data.  The default value is 1000 samples. }
  \item{ncdf.filename}{ A string that will be used as the name of the netCDF file.  This is to be used later when opening and retrieving data generated from this function. }
  \item{silent}{ Logical value.  If \code{FALSE}, the function returns a table of genotype counts generated.  The default is \code{TRUE}; no data will be returned in this case.}  
}
\details{
  The resulting netCDF file will have the following characteristics:
  
Dimensions:

'snp': n.snps*n.chromosomes length

'sample': n.samples length 


Variables:

'sampleID': sample dimension, values 1-n.samples

'position': snp dimension, values [1,2,...,n.chromosomes] n.snps times

'chromosome': snp dimension, values [1,1,...]n.snps times, [2,2,...]n.snps times, ..., [n.chromosomes,n.chromosomes,...]n.snps times

'genotype': 2-dimensional snp x sample, values 0, 1, 2 chosen from allele frequencies that were generated from a uniform distribution on (0,1).  
The missing rate is 0.05 (constant across all SNPs) and is denoted by -1. 
}
\value{
This function returns a table of genotype calls if the silent variable is set to \code{FALSE}, where 0 indicates an AA genotype, 1 is AB, 2 is BB and -1 corresponds to a missing genotype call.

A netCDF file is created from this function and written to disk.  This file (and data) can be accessed later by using the command \code{open.ncdf(ncdf.filename)}.
}
\author{ Caitlin McHugh }
\seealso{\code{\link{ncdf}}, \code{\link{missingGenotypeBySnpSex}}, \code{\link{missingGenotypeByScanChrom}}, \code{\link{simulateIntensityMatrix}} }
\examples{
filenm <- tempfile()

simulateGenotypeMatrix(ncdf.filename=filenm )

file <- NcdfGenotypeReader(filenm)
file #notice the dimensions and variables listed

genot <- getGenotype(file)
table(genot) #can see the number of missing calls

chrom <- getChromosome(file)
unique(chrom) #there are indeed 10 chromosomes, as specified in the function call

close(file)
unlink(filenm)
}
\keyword{datagen}
