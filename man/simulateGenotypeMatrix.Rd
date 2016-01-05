\name{simulateGenotypeMatrix}
\alias{simulateGenotypeMatrix}
\alias{simulateIntensityMatrix}
\title{ Simulate Genotype or Intensity Matrix & Load into GDS/NetCDF File}
\description{
These functions create a simulated genotype or intensity file for test and examples. 
}
\usage{
simulateGenotypeMatrix(n.snps=10, n.chromosomes=10,
                       n.samples=1000, filename,
                       file.type=c("gds", "ncdf"), silent=TRUE)
}
\arguments{
  \item{n.snps}{ An integer corresponding to the number of SNPs per chromosome, the default value is 10.  For this function, the number of SNPs is assumed to be the same for every chromosome. }
  \item{n.chromosomes}{ An integer value describing the total number of chromosomes with default value 10. }
  \item{n.samples}{ An integer representing the number of samples for our data.  The default value is 1000 samples. }
  \item{filename}{ A string that will be used as the name of the file.  This is to be used later when opening and retrieving data generated from this function. }
  \item{file.type}{The type of file to create ("gds" or "ncdf")}
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

OR

'quality': 2-dimensional snp x sample, values between 0 and 1 chosen randomly from a uniform distribution.  There is one quality value per snp, so this value is constant across all samples.

'X': 2-dimensional snp x sample, value of X intensity taken from a normal distribution.  The mean of the distribution for each SNP is based upon the sample genotype.  Mean is 0,2 if sample is homozygous, 1 if heterozygous.

'Y': 2-dimensional snp x sample, value of Y intensity also chosen from a normal distribution, where the mean is chosen according to the mean of X so that sum of means = 2.

}
\value{
\code{simulateGenotypeMatrix} returns a table of genotype calls if the silent variable is set to \code{FALSE}, where 2 indicates an AA genotype, 1 is AB, 0 is BB and -1 corresponds to a missing genotype call.

\code{simulateIntensityMatrix} returns a list if the silent variable is set to FALSE,
which includes:
\item{het}{Heterozygosity table}
\item{nmiss}{Number of missing values}

A file is created and written to disk. 
}
\author{ Caitlin McHugh }
\seealso{\code{\link{GdsGenotypeReader}}, \code{\link{GdsIntensityReader}},
\code{\link{NcdfGenotypeReader}}, \code{\link{NcdfIntensityReader}}}
\examples{
filenm <- tempfile()

simulateGenotypeMatrix(filename=filenm )

file <- GdsGenotypeReader(filenm)
file #notice the dimensions and variables listed

genot <- getGenotype(file)
table(genot) #can see the number of missing calls

chrom <- getChromosome(file)
unique(chrom) #there are indeed 10 chromosomes, as specified in the function call

close(file)

simulateIntensityMatrix(filename=filenm, silent=FALSE )

file <- GdsIntensityReader(filenm)
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
