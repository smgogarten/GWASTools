\name{GenotypeIterator-class}
\docType{class}

\alias{GenotypeIterator}
\alias{GenotypeIterator-class}
\alias{GenotypeBlockIterator}
\alias{GenotypeBlockIterator-class}
\alias{getSnpID,GenotypeIterator-method}
\alias{getChromosome,GenotypeIterator-method}
\alias{getPosition,GenotypeIterator-method}
\alias{getAlleleA,GenotypeIterator-method}
\alias{getAlleleB,GenotypeIterator-method}
\alias{getSnpVariable,GenotypeIterator-method}
\alias{getGenotypeSelection,GenotypeIterator-method}
\alias{iterateFilter}
\alias{iterateFilter,GenotypeIterator-method}
\alias{lastFilter}
\alias{lastFilter,GenotypeIterator-method}
\alias{lastFilter<-}
\alias{lastFilter<-,GenotypeIterator,numeric-method}
\alias{snpFilter}
\alias{snpFilter,GenotypeIterator-method}
\alias{currentFilter}
\alias{currentFilter,GenotypeIterator-method}
\alias{resetIterator}
\alias{resetIterator,GenotypeIterator-method}

\title{Class GenotypeIterator}

\description{
	Extends \code{\link{GenotypeData}} to provide iterators over SNPs.
}

\details{
	Iterator classes allow for iterating over blocks of SNPs.

	For \code{GenotypeBlockIterator}, each call to \code{iterateFilter} will select the next unit of \code{snpBlock} SNPs.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenotypeBlockIterator(genoData, snpBlock=10000, snpInclude=NULL}: Returns a \code{GenotypeBlockIterator} object with the filter set to the first block.

      \code{genoData} is a \code{\link{GenotypeData}} object.

      \code{snpBlock} is an integer specifying the number of SNPs in an iteration block.

      \code{snpInclude} is a vector of snpIDs to include in the filter. If \code{NULL} (default), all SNPs are included.
    }
  }
}

\section{Accessors}{
  \describe{
    \item{}{\code{iterateFilter(x)}: Advance the filter to the next block, range, or set of ranges. Returns \code{TRUE} while there are still variants left to be read, \code{FALSE} if the end of iteration is reached.
    }
    \item{}{\code{lastFilter(x)}, \code{lastFilter(x)<- value}: Get or set the last filter index from the previous call to \code{iterateFilter}.
    }
    \item{}{\code{snpFilter(x)}: Get the list of SNP indices.
    }
    \item{}{\code{currentFilter(x)}: Get the indices of SNPs selected in the current iteration.
    }
    \item{}{\code{resetIterator(x)}: Set the filter to the first block (the same SNPs that are selected when the iterator object is created).
    }
    \item{}{
      \code{getSnpID(object, ...)}: A unique integer vector of snp
      IDs. 
    }
    \item{}{
      \code{getChromosome(object, ...)}: A vector of
      chromosomes. 
    }
    \item{}{
      \code{getPosition(object, ...)}: An integer vector of base pair
      positions.
    }
    \item{}{
      \code{getAlleleA(object, ...)}: A character vector of A alleles.
    }
    \item{}{
      \code{getAlleleB(object, ...)}: A character vector of B alleles.
    }
    \item{}{
      \code{getGenotypeSelection(object, ...)}:
      Extracts genotype values (number of A alleles).
    }
    \item{}{
      \code{getSnpVariable(object, varname, ...)}: Returns the snp
      annotation variable \code{varname}.
    }
    
    Additional arguments for these accessors are passed to the corresponding methods for \code{\link{GenotypeData}} objects.

    Other methods are inherited from \code{\link{GenotypeData}}. 
  }
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{GenotypeData}}
}

\examples{
library(GWASdata)
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)
genoData <- GenotypeData(gds)

iterator <- GenotypeBlockIterator(genoData, snpBlock=100)
id <- getSnpID(iterator)
length(id)
head(id)
geno <- getGenotypeSelection(iterator)
dim(geno)
geno <- getGenotypeSelection(iterator, scan=1:10)
dim(geno)
head(geno)

iterateFilter(iterator)
id <- getSnpID(iterator)
length(id)
head(id)
geno <- getGenotypeSelection(iterator, scan=1:10)
dim(geno)
head(geno)

close(iterator)
}

\keyword{methods}
\keyword{classes}

