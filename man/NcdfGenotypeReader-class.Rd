\name{NcdfGenotypeReader}
\docType{class}

\alias{NcdfGenotypeReader-class}
\alias{NcdfGenotypeReader}
\alias{getVariable,NcdfGenotypeReader-method}
\alias{getSnpID,NcdfGenotypeReader-method}
\alias{getChromosome,NcdfGenotypeReader-method}
\alias{getPosition,NcdfGenotypeReader-method}
\alias{getScanID,NcdfGenotypeReader-method}
\alias{getGenotype,NcdfGenotypeReader-method}
\alias{nsnp,NcdfGenotypeReader-method}
\alias{nscan,NcdfGenotypeReader-method}
\alias{XchromCode,NcdfGenotypeReader-method}
\alias{XYchromCode,NcdfGenotypeReader-method}
\alias{YchromCode,NcdfGenotypeReader-method}
\alias{MchromCode,NcdfGenotypeReader-method}

\title{Class NcdfGenotypeReader}

\description{
  The NcdfGenotypeReader class is an extension of the NcdfReader class
  specific to reading genotype data stored in NetCDF files.
}

\section{Extends}{
  \code{\link{NcdfReader}}
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{NcdfGenotypeReader(filename)}:
      
      \code{filename} must be the path to a NetCDF file.  The NetCDF
      file must contain the following variables:
      \itemize{
	\item 'snp': a coordinate variable with a unique
	integer vector of snp ids
	\item 'chromosome':  integer
	chromosome codes of dimension 'snp'
	\item 'position': integer position values of dimension 'snp'
	\item 'sampleID': a unique integer vector of scan ids with
	dimension 'sample'
	\item 'genotype':
	a matrix of bytes with dimensions ('snp','sample').  The byte values
	must be the number of A alleles : 2=AA, 1=AB, 0=BB.
      }
      Default values for chromosome codes are 1-22, 23=X, 24=XY, 25=Y,
      26=M.  The defaults may be changed with the arguments
      \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
      \code{MchromCode}.
      
      The \code{NcdfGenotypeReader} constructor creates and returns a
      NcdfGenotypeReader instance pointing to this file.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a NcdfGenotypeReader object.
  \code{snp} and
      \code{scan} indicate which elements to return along the snp and
      scan dimensions.  They must be integer vectors of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that the entire dimension should be read.  If snp
      and/or is scan omitted, the entire variable is read.
      
      See \code{\link{NcdfReader}} for additional methods.

  \describe{
    \item{}{
      \code{nsnp(object)}: The number of SNPs in the NetCDF file.
    }
    \item{}{
      \code{nscan(object)}: The number of scans in the NetCDF file.
    }
    \item{}{
      \code{getSnpID(object, index)}: A unique integer vector of snp
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{getChromosome(object, index, char=FALSE)}: A vector of
      chromosomes.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.
      If \code{char=FALSE} (default), returns an integer vector.
      If \code{char=TRUE}, returns a character vector with elements in
      (1:22,X,XY,Y,M,U).  "U" stands for "Unknown" and is the value
      given to any chromosome code not falling in the other categories.
    }
    \item{}{
      \code{getPosition(object, index)}: An integer vector of base pair
      positions.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{getScanID(object, index)}: A unique integer vector of scan
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{getGenotype(object, snp, scan)}: Extracts genotype values
      (number of A alleles).  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{getVariable(object, varname, snp, scan)}: Extracts the
      contents of the variable \code{varname}.  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.  If the variable is not found in the NetCDF
      file, returns \code{NULL}.
    }
     \item{}{
       \code{XchromCode(object)}: Returns the integer code for the X
       chromosome.
     }
     \item{}{
       \code{XYchromCode(object)}: Returns the integer code for the
       pseudoautosomal region.
     }
     \item{}{
       \code{YchromCode(object)}: Returns the integer code for the Y
       chromosome.
     }
     \item{}{
       \code{MchromCode(object)}: Returns the integer code for 
       mitochondrial SNPs.
     }
  }
}

\author{Stephanie Gogarten}

\seealso{\code{\link{NcdfReader}},
  \code{\link{NcdfIntensityReader}},
  \code{\link{GenotypeData}}, \code{\link{IntensityData}}
}

\examples{
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(file)

# dimensions
nsnp(nc)
nscan(nc)

# get snpID and chromosome
snpID <- getSnpID(nc)
chrom <- getChromosome(nc)

# get positions only for chromosome 22
pos22 <- getPosition(nc, index=(chrom == 22))

# get all snps for first scan
geno <- getGenotype(nc, snp=c(1,-1), scan=c(1,1))

# starting at snp 100, get 10 snps for the first 5 scans
geno <- getGenotype(nc, snp=c(100,10), scan=c(1,5))

close(nc)
}

\keyword{methods}
\keyword{classes}

