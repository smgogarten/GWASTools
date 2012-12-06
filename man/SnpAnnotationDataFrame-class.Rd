\name{SnpAnnotationDataFrame}
\docType{class}

\alias{SnpAnnotationDataFrame-class}
\alias{SnpAnnotationDataFrame}
\alias{hasVariable,SnpAnnotationDataFrame-method}
\alias{getVariable,SnpAnnotationDataFrame-method}
\alias{getSnpID,SnpAnnotationDataFrame-method}
\alias{getChromosome,SnpAnnotationDataFrame-method}
\alias{getPosition,SnpAnnotationDataFrame-method}
\alias{getAlleleA,SnpAnnotationDataFrame-method}
\alias{getAlleleB,SnpAnnotationDataFrame-method}
\alias{getVariableNames,SnpAnnotationDataFrame-method}
\alias{getAnnotation,SnpAnnotationDataFrame-method}
\alias{getMetadata,SnpAnnotationDataFrame-method}
\alias{autosomeCode,SnpAnnotationDataFrame-method}
\alias{XchromCode,SnpAnnotationDataFrame-method}
\alias{XYchromCode,SnpAnnotationDataFrame-method}
\alias{YchromCode,SnpAnnotationDataFrame-method}
\alias{MchromCode,SnpAnnotationDataFrame-method}

\title{Class SnpAnotationDataFrame}

\description{
  The SnpAnnotationDataFrame class stores annotation data associated with
  SNPs, as well as metadata describing each column.  It extends the
  \code{\link{AnnotatedDataFrame}} class.
}

\section{Extends}{
  \code{\link{AnnotatedDataFrame}}
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{SnpAnnotationDataFrame(data, metadata)}:

      \code{data} must be a data.frame containing the SNP annotation.
      It must contain at least the following columns:
      \itemize{
	\item "snpID": integer vector containing unique SNP ids.
	\item "chromosome": integer vector containing chromosome codes.
	\item "position": integer vector containing position (in base
	pairs) on the chromosome.
      }
      Default values for chromosome codes are 1-22=autosome, 23=X, 24=XY, 25=Y,
      26=M.  The defaults may be changed with the arguments \code{autosomeCode},
      \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
      \code{MchromCode}.

      \code{metadata} is an optional data.frame containing a description
      for each column in \code{data}.  It should contain a column
      "labelDescription", with \code{row.names(metadata) == names(data)}.

      The \code{SnpAnnotationDataFrame} constructor creates and returns
      a SnpAnnotationDataFrame instance.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a SnpAnnotationDataFrame
  object.

  \describe{
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
      \code{getVariable(object, varname, index)}: A vector of the
      column \code{varname}.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.
      If \code{varname} is itself a vector, returns a data.frame.
      Returns \code{NULL} if
      \code{varname} is not found in \code{object}.
    }
    \item{}{
      \code{hasVariable(object, varname)}: Returns \code{TRUE} if
      \code{varname} is a column in \code{object}, \code{FALSE} if not.
    }
    \item{}{
      \code{getVariableNames(object)}: Returns a character vector with
      the names of all columns in \code{object}.
    }
    \item{}{
      \code{getAnnotation(object)}: Returns all annotation variables
      as a data frame.
    }
    \item{}{
      \code{getMetadata(object)}: Returns metadata describing the
      annotation variables as a data frame.
    }
    
    Inherited methods from
    \code{\link{AnnotatedDataFrame}}:
    \item{}{
      \code{varLabels(object)}: Returns a character vector with
      the names of all columns in \code{object}.
    }
    \item{}{
      \code{pData(object)}: Returns all annotation variables
      as a data frame, or sets the annotation variables with
      \code{pData(object) <- df}.
    }
    \item{}{
      \code{varMetadata(object)}: Returns metadata describing the
      annotation variables as a data frame, or sets the metadata with
      \code{varMetadata(object) <- df}.
    }
    \item{}{
      The operators \code{[}, \code{$}, and \code{[[} work just as they do in
      standard data frames, for both retrieval and assignment.
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
  \code{\link{AnnotatedDataFrame}}, \code{\link{ScanAnnotationDataFrame}},
  \code{\link{GenotypeData}}, \code{\link{IntensityData}}
}

\examples{
library(GWASdata)
data(affy_snp_annot)
snpAnnot <- SnpAnnotationDataFrame(affy_snp_annot)

# list columns
varLabels(snpAnnot)

# add metadata
meta <- varMetadata(snpAnnot)
meta["snpID", "labelDescription"] <- "unique integer ID"
varMetadata(snpAnnot) <- meta

# get snpID and chromosome
snpID <- getSnpID(snpAnnot)
chrom <- getChromosome(snpAnnot)

# get positions only for chromosome 22
pos22 <- getPosition(snpAnnot, index=(chrom == 22))

# get rsID
if (hasVariable(snpAnnot, "rsID")) rsID <- getVariable(snpAnnot, "rsID")

# display data
head(pData(snpAnnot))

# standard operators
snpID <- snpAnnot$snpID
chrom <- snpAnnot[["chromosome"]]
subset <- snpAnnot[1:10, 1:5]
snpAnnot$newVar <- rep(1, nrow(snpAnnot))

# replace data
df <- pData(snpAnnot)
pData(snpAnnot) <- df

# PLINK chromosome coding
snpID <- 1:10
chrom <- c(rep(1L,5), 23:27)
pos <- 101:110
df <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
snpAnnot <- SnpAnnotationDataFrame(df, YchromCode=24L, XYchromCode=25L)
getChromosome(snpAnnot, char=TRUE)
}

\keyword{methods}
\keyword{classes}
    