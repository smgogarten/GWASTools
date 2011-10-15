\name{ScanAnnotationDataFrame}
\docType{class}

\alias{ScanAnnotationDataFrame-class}
\alias{ScanAnnotationDataFrame}
\alias{hasVariable,ScanAnnotationDataFrame-method}
\alias{getVariable,ScanAnnotationDataFrame-method}
\alias{getScanID,ScanAnnotationDataFrame-method}
\alias{hasSex,ScanAnnotationDataFrame-method}
\alias{getSex,ScanAnnotationDataFrame-method}
\alias{getVariableNames,ScanAnnotationDataFrame-method}
\alias{getAnnotation,ScanAnnotationDataFrame-method}
\alias{getMetadata,ScanAnnotationDataFrame-method}

\title{Class ScanAnotationDataFrame}

\description{
  The ScanAnnotationDataFrame class stores annotation data associated with
  subjects in a genotyping study, where there may be multiple scans per
  subject, as well as metadata describing each column.  It extends the
  \code{\link{AnnotatedDataFrame}} class.
}

\section{Extends}{
  \code{\link{AnnotatedDataFrame}}
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{ScanAnnotationDataFrame(data, metadata)}:

      \code{data} must be a data.frame containing the scan annotation.
      It must contain at least the following column:
      \itemize{
	\item "scanID": integer vector containing unique scan ids.
      }
      If a column representing sex is present, it must have the
      following format:
      \itemize{
	\item "sex": character vector with values 'M' or 'F'.
      }

      \code{metadata} is an optional data.frame containing a description
      for each column in \code{data}.  It should contain a column
      "labelDescription", with \code{row.names(metadata) == names(data)}.
      
      The \code{ScanAnnotationDataFrame} constructor creates and returns
      a ScanAnnotationDataFrame instance.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a ScanAnnotationDataFrame
  object.

  \describe{
    \item{}{
      \code{getScanID(object, index)}: A unique integer vector of scan
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{getSex(object, index)}: A character vector of sex, with values 'M'
      or 'F'.    The optional \code{index} is a logical or
      integer vector specifying elements to extract.  
    }
    \item{}{
      \code{hasSex(object)}: Returns \code{TRUE} if the column 'sex' is present in
      \code{object}. 
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
      The operators \code{$} and \code{[} work just as they do in
      standard data frames, for both retrieval and assignment.
    }
  }
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{AnnotatedDataFrame}}, \code{\link{SnpAnnotationDataFrame}},
  \code{\link{GenotypeData}}, \code{\link{IntensityData}}
}

\examples{
library(GWASdata)
data(affy_scan_annot)
scanAnnot <- ScanAnnotationDataFrame(affy_scan_annot)

scanID <- getScanID(scanAnnot)
sex <- getSex(scanAnnot)
if (hasVariable(scanAnnot, "plate")) plate <- getVariable(scanAnnot, "plate")
subjectID <- getVariable(scanAnnot, "subjectID", index=(sex == "M"))

# list columns
varLabels(scanAnnot)

# add metadata
meta <- varMetadata(scanAnnot)
meta["scanID", "labelDescription"] <- "unique integer ID"
varMetadata(scanAnnot) <- meta

# display data
head(pData(scanAnnot))

# standard operators
scanID <- scanAnnot$scanID
sex <- scanAnnot[["sex"]]
subset <- scanAnnot[1:10, 1:5]
scanAnnot$newVar <- rep(1, nrow(scanAnnot))

# replace data
df <- pData(scanAnnot)
pData(scanAnnot) <- df
}

\keyword{methods}
\keyword{classes}
    