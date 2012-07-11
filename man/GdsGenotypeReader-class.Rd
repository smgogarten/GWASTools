\name{GdsGenotypeReader}
\docType{class}

\alias{GdsGenotypeReader-class}
\alias{GdsGenotypeReader}
\alias{getVariable,GdsGenotypeReader-method}
\alias{getSnpID,GdsGenotypeReader-method}
\alias{getChromosome,GdsGenotypeReader-method}
\alias{getPosition,GdsGenotypeReader-method}
\alias{getScanID,GdsGenotypeReader-method}
\alias{getGenotype,GdsGenotypeReader-method}
\alias{nsnp,GdsGenotypeReader-method}
\alias{nscan,GdsGenotypeReader-method}
\alias{XchromCode,GdsGenotypeReader-method}
\alias{XYchromCode,GdsGenotypeReader-method}
\alias{YchromCode,GdsGenotypeReader-method}
\alias{MchromCode,GdsGenotypeReader-method}

\title{Class GdsGenotypeReader}

\description{
  The GdsGenotypeReader class is an extension of the GdsReader class
  specific to reading genotype data stored in GDS files.
}

\section{Extends}{
  \code{\link{GdsReader}}
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GdsGenotypeReader(filename)}:
      
      \code{filename} must be the path to a GDS file.  The GDS
      file must contain the following variables:
      \itemize{
	\item 'snp.id': a unique integer vector of snp ids
	\item 'snp.chromosome':  integer chromosome codes
	\item 'snp.position': integer position values
	\item 'sample.id': a unique integer vector of scan ids
	\item 'genotype':
	a matrix of bytes with dimensions ('snp','sample').  The byte values
	must be the number of A alleles : 2=AA, 1=AB, 0=BB.
      }
      Default values for chromosome codes are 1-22, 23=X, 24=XY, 25=Y,
      26=M.  The defaults may be changed with the arguments
      \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
      \code{MchromCode}.
      
      The \code{GdsGenotypeReader} constructor creates and returns a
      GdsGenotypeReader instance pointing to this file.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a GdsGenotypeReader object.
  \code{snp} and
      \code{scan} indicate which elements to return along the snp and
      scan dimensions.  They must be integer vectors of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that the entire dimension should be read.  If snp
      and/or is scan omitted, the entire variable is read.
      
      See \code{\link{GdsReader}} for additional methods.

  \describe{
    \item{}{
      \code{nsnp(object)}: The number of SNPs in the GDS file.
    }
    \item{}{
      \code{nscan(object)}: The number of scans in the GDS file.
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
      represented as \code{NA}.  If the variable is not found in the GDS
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

\seealso{\code{\link{GdsReader}},
  \code{\link{GenotypeData}}
}

\examples{
library(SNPRelate)
gds <- GdsGenotypeReader(snpgdsExampleFileName())

# dimensions
nsnp(gds)
nscan(gds)

# get snpID and chromosome
snpID <- getSnpID(gds)
chrom <- getChromosome(gds)

# get positions only for chromosome 22
pos22 <- getPosition(gds, index=(chrom == 22))

# get all snps for first scan
geno <- getGenotype(gds, snp=c(1,-1), scan=c(1,1))

# starting at snp 100, get 10 snps for the first 5 scans
geno <- getGenotype(gds, snp=c(100,10), scan=c(1,5))

close(gds)
}

\keyword{methods}
\keyword{classes}

