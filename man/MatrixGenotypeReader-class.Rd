\name{MatrixGenotypeReader}
\docType{class}

\alias{MatrixGenotypeReader-class}
\alias{MatrixGenotypeReader}
\alias{getSnpID,MatrixGenotypeReader-method}
\alias{getChromosome,MatrixGenotypeReader-method}
\alias{getPosition,MatrixGenotypeReader-method}
\alias{getScanID,MatrixGenotypeReader-method}
\alias{getGenotype,MatrixGenotypeReader-method}
\alias{nsnp,MatrixGenotypeReader-method}
\alias{nscan,MatrixGenotypeReader-method}
\alias{show,MatrixGenotypeReader-method}
\alias{autosomeCode,MatrixGenotypeReader-method}
\alias{XchromCode,MatrixGenotypeReader-method}
\alias{XYchromCode,MatrixGenotypeReader-method}
\alias{YchromCode,MatrixGenotypeReader-method}
\alias{MchromCode,MatrixGenotypeReader-method}

\title{Class MatrixGenotypeReader}

\description{
  The MatrixGenotypeReader class stores a matrix of genotypes as well as
  SNP and scan IDs, chromosome, and position.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{MatrixGenotypeReader(genotype=genotype, snpID=snpID, chromosome=chromosome, position=position, scanID=scanID)}:
      
      \code{genotype} must be a matrix with dimensions ('snp','scan') containing the number of A alleles : 2=AA, 1=AB, 0=BB.
      
      \code{snp} must be a unique integer vector of SNP ids.
      
      \code{chromosome} must be an integer vector of chromosomes.
      Default values for chromosome codes are 1-22=autosome, 23=X, 24=XY, 25=Y,
      26=M.  The defaults may be changed with the arguments \code{autosomeCode},
      \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
      \code{MchromCode}.
      
      \code{position} must be an integer vector of base positions
      
      \code{scanID} must be a unique integer vector of scan ids .
      }
      
      The \code{MatrixGenotypeReader} constructor creates and returns a
      MatrixGenotypeReader instance.
    }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a MatrixGenotypeReader object.

  \describe{
    \item{}{
      \code{nsnp(object)}: The number of SNPs.
    }
    \item{}{
      \code{nscan(object)}: The number of scans.
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
      \code{getGenotype(object, snp=c(1,-1), scan=c(1,-1))}: Extracts genotype values
      (number of A alleles).  
      \code{snp} and 
      \code{scan} indicate which elements to return along the snp and
      scan dimensions.  They must be integer vectors of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that the entire dimension should be read.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
     \item{}{
       \code{autosomeCode(object)}: Returns the integer codes for the
       autosomes.
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

\seealso{
  \code{\link{NcdfGenotypeReader}},
  \code{\link{GenotypeData}}
}

\examples{
snpID <- 1:100
chrom <- rep(1:20, each=5)
pos <- 1001:1100
scanID <- 1:20
geno <- matrix(sample(c(0,1,2,NA), 2000, replace=TRUE), nrow=100, ncol=20)

mgr <- MatrixGenotypeReader(genotype=geno, snpID=snpID,
  chromosome=chrom, position=pos, scanID=scanID)

# dimensions
nsnp(mgr)
nscan(mgr)

# get snpID and chromosome
snpID <- getSnpID(mgr)
chrom <- getChromosome(mgr)

# get positions only for chromosome 10
pos10 <- getPosition(mgr, index=(chrom == 10))

# get all snps for first scan
geno <- getGenotype(mgr, snp=c(1,-1), scan=c(1,1))

# starting at snp 50, get 10 snps for the first 5 scans
geno <- getGenotype(mgr, snp=c(50,10), scan=c(1,5))
}

\keyword{methods}
\keyword{classes}

