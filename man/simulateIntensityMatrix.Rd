\name{simulateIntensityMatrix}
\alias{simulateIntensityMatrix}
\title{Simulate Intensity Matrix & Load into NetCDF File}
\description{
This function creates a netCDF file with dimensions 'snp' and 'sample'
and variables 'sampleID', 'position', 'chromosome', 'quality', 'X', and 'Y'.  These variables hold simulated data as explained below.  Mainly, this function is intended to be used in examples involving matrices holding quantitative data.
}
\usage{
simulateIntensityMatrix(n.snps=10, n.chromosomes=10,
                        n.samples=1000, ncdf.filename,
                        silent=TRUE)
}
\arguments{
  \item{n.snps}{ An integer corresponding to the number of SNPs per chromosome, the default value is 10.  For this function, the number of SNPs is assumed to be the same for every chromosome. }
  \item{n.chromosomes}{ An integer value describing the total number of chromosomes with default value 10. }
  \item{n.samples}{ An integer representing the number of samples for our data.  The default value is 1000 samples. }
  \item{ncdf.filename}{ A string that will be used as the name of the netCDF file.  This is to be used later when opening and retrieving data generated from this function. }
  \item{silent}{ Logical value.  If \code{FALSE}, the function returns a
    list of heterozygosity and missing values.  The default is \code{TRUE}; no data will be returned in this case. }
}
\details{
  The resulting netCDF file will have the following characteristics:
  
Dimensions:

'snp': n.snps*n.chromosomes length

'sample': n.samples length


Variables:

'sampleID': sample dimension, values 1-n.samples

'position': snp dimension, values [1,2,...,n.chromosomes] n.snps times

'chromosome': snp dimension, values[1,1,...]n.snps times, [2,2,...]n.snps times, ... , [n.chromosomes,n.chromosomes,...]n.snps times

'quality': 2-dimensional snp x sample, values between 0 and 1 chosen randomly from a uniform distribution.  There is one quality value per snp, so this value is constant across all samples.

'X': 2-dimensional snp x sample, value of X intensity taken from a normal distribution.  The mean of the distribution for each SNP is based upon the sample genotype.  Mean is 0,2 if sample is homozygous, 1 if heterozygous.

'Y': 2-dimensional snp x sample, value of Y intensity also chosen from a normal distribution, where the mean is chosen according to the mean of X so that sum of means = 2.

}
\value{
This function returns a list if the silent variable is set to FALSE,
which includes:
\item{het}{Heterozygosity table}
\item{nmiss}{Number of missing values}

A netCDF file is created from this function and written to disk.  This file (and data) can be accessed later by using the command 'open.ncdf(ncdf.filename)'.
}
\author{ Caitlin McHugh }
\seealso{\code{\link{ncdf}}, \code{\link{meanIntensityByScanChrom}}, \code{\link{simulateGenotypeMatrix}} }
\examples{
filenm <- tempfile()

simulateIntensityMatrix(ncdf.filename=filenm, silent=FALSE )

file <- NcdfIntensityReader(filenm)
file #notice the dimensions and variables listed

xint <- getX(file)
yint <- getY(file)
print("Number missing is: "); sum(is.na(xint))

chrom <- getChromosome(file)
unique(chrom) #there are indeed 10 chromosomes, as specified in the function call

close(file)
unlink(filenm)

}
\keyword{datagen}
