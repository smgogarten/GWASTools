\name{GdsGenotypeReader}
\docType{class}

\alias{GdsGenotypeReader-class}
\alias{GdsGenotypeReader}
\alias{getVariable,GdsGenotypeReader-method}
\alias{getSnpID,GdsGenotypeReader-method}
\alias{getChromosome,GdsGenotypeReader-method}
\alias{getPosition,GdsGenotypeReader-method}
\alias{getAlleleA,GdsGenotypeReader-method}
\alias{getAlleleB,GdsGenotypeReader-method}
\alias{getScanID,GdsGenotypeReader-method}
\alias{getGenotype,GdsGenotypeReader-method}
\alias{getGenotypeSelection,GdsGenotypeReader-method}
\alias{nsnp,GdsGenotypeReader-method}
\alias{nscan,GdsGenotypeReader-method}
\alias{autosomeCode,GdsGenotypeReader-method}
\alias{XchromCode,GdsGenotypeReader-method}
\alias{XYchromCode,GdsGenotypeReader-method}
\alias{YchromCode,GdsGenotypeReader-method}
\alias{MchromCode,GdsGenotypeReader-method}

\title{Class GdsGenotypeReader}

\description{
  The GdsGenotypeReader class is an extension of the GdsReader class
  specific to reading genotype data stored in GDS files. GDS files with both \code{snp x scan} and \code{scan x snp} dimensions are supported.
}

\section{Extends}{
  \code{\link{GdsReader}}
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GdsGenotypeReader(filename, genotypeDim)}:
      
      \code{filename} must be the path to a GDS file or a gds object.  
      The GDS file must contain the following variables:
      \itemize{
	\item 'snp.id': a unique integer vector of snp ids
	\item 'snp.chromosome':  integer chromosome codes
	\item 'snp.position': integer position values
	\item 'sample.id': a unique integer vector of scan ids
	\item 'genotype':
	a matrix of bytes with dimensions ('snp','sample').  The byte values
	must be the number of A alleles : 2=AA, 1=AB, 0=BB.

      }
        
      The optional variable "snp.allele" stores the A and B alleles in a
      character vector with format "A/B".
      
      Default values for chromosome codes are 1-22=autosome, 23=X, 24=XY, 25=Y,
      26=M.  The defaults may be changed with the arguments \code{autosomeCode},
      \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
      \code{MchromCode}.

      The constructor automatically detects whether the GDS file is in \code{snp x scan} or \code{scan x snp} order using the dimensions of \code{snp.id} and \code{sample.id}. In the case of GDS files with equal SNP and scan dimensions, \code{genotypeDim} is a required input to the function and can take values \code{"snp,scan"} or \code{"scan,snp"}.

      The \code{GdsGenotypeReader} constructor creates and returns a
      GdsGenotypeReader instance pointing to this file.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a GdsGenotypeReader object.     
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
      \code{getAlleleA(object, index)}: A character vector of A alleles.
      The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{getAlleleB(object, index)}: A character vector of B alleles.
      The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{getScanID(object, index)}: A unique integer vector of scan
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{getGenotype(object, snp=c(1,-1), scan=c(1,-1), drop=TRUE, transpose=FALSE, ...)}: Extracts genotype values
      (number of A alleles).  
  \code{snp} and
      \code{scan} indicate which elements to return along the snp and
      scan dimensions.  They must be integer vectors of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that the entire dimension should be read.
      If \code{drop=TRUE} the result is coerced to the lowest possible dimension.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values and the value of \code{drop}.  Missing values are
      represented as \code{NA}. Genotypes are returned in SNP x scan order if transpose=FALSE, otherwise they are returned in scan x SNP order.
    }
    \item{}{
      \code{getGenotypeSelection(object, snp=NULL, scan=NULL, drop=TRUE, transpose=FALSE, ...)}: Extracts genotype values
      (number of A alleles).  
      \code{snp} and
      \code{scan} may be integer or logical vectors indicating which elements to return along the snp and
      scan dimensions.  Unlike \code{getGenotype}, the values requested need not be in contiguous blocks.
      Other arguments are identical to \code{getGenotype}.
    }
    \item{}{
      \code{getVariable(object, varname, index, drop=TRUE, ...)}: Extracts the
      contents of the variable \code{varname}.    
      The optional \code{index} is a logical or
      integer vector (if \code{varname} is 1D) or list (if \code{varname} is 2D or higher) specifying elements to extract.
      If \code{drop=TRUE} the result is coerced to the lowest possible dimension.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values and the value of \code{drop}.  Missing values are
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
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)

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

# get snps 1-10, 25-30 for scans 3,5,7
snp.index <- c(1:10, 25:30)
scan.index <- c(3,5,7)
geno <- getGenotypeSelection(gds, snp=snp.index, scan=scan.index)
dimnames(geno) <- list(getSnpID(gds, index=snp.index),
                       getScanID(gds, index=scan.index))

close(gds)
}

\keyword{methods}
\keyword{classes}

